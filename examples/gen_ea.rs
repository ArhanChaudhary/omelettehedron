//! Table design study: generate the edge x A-centre table (`optimal::ea`) and save it.
//!
//! usage: cargo run --release --example gen_ea -- out.prun

use omelettehedron::optimal::ea::{build_ea, EaIndex};
use std::sync::atomic::AtomicBool;

fn main() -> std::io::Result<()> {
    let out = std::env::args().nth(1).unwrap_or_else(|| "opt-ea.prun".into());
    let t0 = std::time::Instant::now();
    let ix = EaIndex::new()?;
    eprintln!(
        "index: {} classes, {} entries ({:.1} GB), {:.1}s",
        ix.rep.len(),
        ix.size(),
        ix.size() as f64 / 4e9,
        t0.elapsed().as_secs_f64()
    );
    let stop = AtomicBool::new(false);
    let t = build_ea(&ix, true, &stop).ok_or_else(|| std::io::Error::other("stopped"))?;
    eprintln!("generated in {:.1}s; histogram {:?}", t0.elapsed().as_secs_f64(), t.hist);
    t.save(std::path::Path::new(&out))?;
    eprintln!("saved {out} at {:.1}s", t0.elapsed().as_secs_f64());
    Ok(())
}
