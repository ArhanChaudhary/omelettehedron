//! Analysis-only pruning statistics (feature `stats`): per remaining budget, which bounds
//! would prune each evaluated child. Thread-local counters, collected with `rayon::broadcast`.

use std::cell::RefCell;

pub const NB: usize = 24;
/// condition bits
pub const AE: usize = 1;
pub const CB: usize = 2;
pub const CBE: usize = 4;
pub const BAE: usize = 8;
pub const BCB: usize = 16;
/// hypothetical: CBE on the class-B position (not used for pruning)
pub const BCBE: usize = 32;
/// hypothetical: forced-route CB bound (distance to the CB projection of
/// S_s = {L in H_A : d2(L) <= s}, s = best - 1 - iteration bound)
pub const FR: usize = 64;
pub const NPAT: usize = 128;

/// Forced-route CB tables for s = 0..=7 (exact distances, 24-fold symmetry-reduced index).
pub struct Route {
    pub ix: crate::cbsym::CbSymIndex,
    pub d: Vec<Vec<u8>>,
}

pub fn route() -> &'static Route {
    static R: std::sync::OnceLock<Route> = std::sync::OnceLock::new();
    R.get_or_init(|| {
        let t0 = std::time::Instant::now();
        let ix = crate::cbsym::CbSymIndex::new();
        let p2 = crate::tables::Phase2::new();
        let attach = crate::tables::CbTable::attach_targets();
        let ne = crate::tables::N_P2_EDGE;
        let cmin: Vec<u8> = (0..crate::coord::N_CORNER)
            .map(|c| p2.dist[c * ne..(c + 1) * ne].iter().copied().min().unwrap_or(99))
            .collect();
        let d: Vec<Vec<u8>> = (0..8)
            .map(|s| {
                let cs: Vec<usize> = (0..crate::coord::N_CORNER)
                    .filter(|&c| cmin[c] as usize <= s)
                    .collect();
                let v = ix.exact_to(&cs, &attach);
                let mut h = [0u64; 16];
                for &x in &v {
                    h[(x as usize).min(15)] += 1;
                }
                let tot: u64 = h.iter().sum();
                let mean = h.iter().enumerate().map(|(i, &c)| i as f64 * c as f64).sum::<f64>() / tot as f64;
                eprintln!("route s={s}: {} target corners, hist {:?} mean {:.3}", cs.len(), h, mean);
                v
            })
            .collect();
        eprintln!("route tables built in {:.1}s", t0.elapsed().as_secs_f64());
        Route { ix, d }
    })
}

#[derive(Clone)]
pub struct Stats {
    /// [budget][pattern] children evaluated with this set of failing conditions
    pub pat: Vec<[u64; NPAT]>,
    /// [budget][slack 0..8] survivors by slack of the phase-1 bound (budget - max(ha,hc,cbe))
    pub slack: Vec<[u64; 8]>,
    /// [budget][slack 0..8] survivors by slack of the class-B bound (lim - max(hba,hbc))
    pub bslack: Vec<[u64; 8]>,
    /// [budget][binding] survivors by which phase-1 bound is the max (0 AE,1 CB,2 CBE, 3 tie)
    pub bind: Vec<[u64; 4]>,
    /// leaves: phase-2 length histogram; [0..16] useful (total < best), [16..32] not
    pub leaf_d2: [u64; 32],
    pub leaf_nonmember: u64,
    /// children evaluated / surviving by allowance = best - 1 - depth (moves left for a
    /// solution shorter than the best); only with a best solution
    pub lim_eval: [u64; NB],
    pub lim_surv: [u64; NB],
    /// [budget] survivors that the forced-route bound would prune
    pub fr_surv: [u64; NB],
}

impl Default for Stats {
    fn default() -> Self {
        Stats {
            pat: vec![[0; NPAT]; NB],
            slack: vec![[0; 8]; NB],
            bslack: vec![[0; 8]; NB],
            bind: vec![[0; 4]; NB],
            leaf_d2: [0; 32],
            leaf_nonmember: 0,
            lim_eval: [0; NB],
            lim_surv: [0; NB],
            fr_surv: [0; NB],
        }
    }
}

thread_local! {
    pub static ST: RefCell<Stats> = RefCell::new(Stats::default());
}

pub fn with<F: FnOnce(&mut Stats)>(f: F) {
    ST.with(|s| f(&mut s.borrow_mut()));
}

/// Collect and reset the counters of every pool thread.
pub fn take() -> Stats {
    let parts: Vec<Stats> = rayon::broadcast(|_| ST.with(|s| std::mem::take(&mut *s.borrow_mut())));
    let mut t = Stats::default();
    for p in parts {
        for b in 0..NB {
            for k in 0..NPAT {
                t.pat[b][k] += p.pat[b][k];
            }
            for k in 0..8 {
                t.slack[b][k] += p.slack[b][k];
                t.bslack[b][k] += p.bslack[b][k];
            }
            for k in 0..4 {
                t.bind[b][k] += p.bind[b][k];
            }
        }
        for k in 0..32 {
            t.leaf_d2[k] += p.leaf_d2[k];
        }
        for k in 0..NB {
            t.fr_surv[k] += p.fr_surv[k];
            t.lim_eval[k] += p.lim_eval[k];
            t.lim_surv[k] += p.lim_surv[k];
        }
        t.leaf_nonmember += p.leaf_nonmember;
    }
    t
}

pub fn print(t: &Stats) {
    eprintln!("stats: budget eval surv | pruned: onlyAE onlyCB onlyCBE onlyB(ae/cb) multi | surv-with-bCBE-prune | slack0 slack1 slack2+ | bind AE/CB/CBE/tie | bslack0 bslack1 bslack2+");
    for b in (0..NB).rev() {
        let row = &t.pat[b];
        let eval: u64 = row.iter().sum();
        if eval == 0 {
            continue;
        }
        let real = |k: usize| k & (AE | CB | CBE | BAE | BCB);
        let surv: u64 = (0..NPAT).filter(|&k| real(k) == 0).map(|k| row[k]).sum();
        let only = |m: usize| -> u64 {
            (0..NPAT)
                .filter(|&k| real(k) != 0 && real(k) & !m == 0)
                .map(|k| row[k])
                .sum()
        };
        let o_ae = only(AE);
        let o_cb = only(CB);
        let o_cbe = only(CBE);
        let o_b = only(BAE | BCB);
        let pruned = eval - surv;
        let multi = pruned - o_ae - o_cb - o_cbe - o_b;
        let sb: u64 = (0..NPAT).filter(|&k| real(k) == 0 && k & BCBE != 0).map(|k| row[k]).sum();
        let s = &t.slack[b];
        let bs = &t.bslack[b];
        let bi = &t.bind[b];
        eprintln!(
            "stats: {:>2} {:>13} {:>13} ({:5.3}) | {:>12} {:>12} {:>12} {:>12} {:>12} | {:>11} | {:>12} {:>12} {:>12} | {} {} {} {} | {} {} {}",
            b, eval, surv, surv as f64 / eval as f64, o_ae, o_cb, o_cbe, o_b, multi, sb,
            s[0], s[1], s[2..].iter().sum::<u64>(), bi[0], bi[1], bi[2], bi[3],
            bs[0], bs[1], bs[2..].iter().sum::<u64>()
        );
    }
    eprintln!("stats: forced-route CB would prune survivors by budget: {:?}", &t.fr_surv);
    eprintln!("stats: by allowance (evaluated, surviving):");
    for a in (0..NB).rev() {
        if t.lim_eval[a] > 0 {
            eprintln!("stats: allow {:>2} {:>14} {:>14}", a, t.lim_eval[a], t.lim_surv[a]);
        }
    }
    eprintln!(
        "stats: leaves useful d2 {:?} useless d2 {:?} nonmember {}",
        &t.leaf_d2[..16],
        &t.leaf_d2[16..],
        t.leaf_nonmember
    );
}
