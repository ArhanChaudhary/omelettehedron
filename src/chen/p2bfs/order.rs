//! Orders of the classes for the search's parts: a part is a range of consecutive classes of
//! the order ([`super::instance::Instance::order`]), and the push of a part reads the layer
//! files of every class with a move into it, so an order whose consecutive classes share their
//! neighbours saves reads.
//!
//! * [`a_move_order`]: the classes joined by the always-allowed moves (U, F, BR, BL) together,
//!   for the block split's space, where those moves never leave a component;
//! * [`grown_order`]: clusters of a given size grown greedily, each next class the one with the
//!   most moves into the cluster so far (a graph-growing partition), for the radius-1 space,
//!   where no move is always allowed;
//! * [`reads`]: the files a layer's step reads for parts of a given size, to compare orders.

use super::instance::NO_EDGE;
use crate::cube::N_MOVES;

/// The neighbours of each class (the classes with a move from or into it, itself included
/// when a move keeps the class), sorted.
pub fn adjacency(edges: &[[(u16, u8); N_MOVES]]) -> Vec<Vec<u32>> {
    let n = edges.len();
    let mut adj: Vec<Vec<u32>> = vec![vec![]; n];
    for (c, row) in edges.iter().enumerate() {
        for &(c2, _) in row {
            if c2 != NO_EDGE {
                adj[c].push(c2 as u32);
                adj[c2 as usize].push(c as u32);
            }
        }
    }
    for a in &mut adj {
        a.sort_unstable();
        a.dedup();
    }
    adj
}

/// The classes ordered by their components under the always-allowed moves (U, F, BR and BL):
/// the components in the order of their smallest class, each in class order.
pub fn a_move_order(edges: &[[(u16, u8); N_MOVES]], pmove: &[u8; N_MOVES]) -> Vec<u32> {
    fn find(root: &mut [usize], mut i: usize) -> usize {
        while root[i] != i {
            root[i] = root[root[i]];
            i = root[i];
        }
        i
    }
    let n = edges.len();
    let mut root: Vec<usize> = (0..n).collect();
    for c in 0..n {
        for m in 0..N_MOVES {
            let (c2, _) = edges[c][m];
            let face = usize::from(pmove[m] / 2);
            if c2 != NO_EDGE && crate::geometry::A_FACES.contains(&face) {
                let (a, b) = (find(&mut root, c), find(&mut root, c2 as usize));
                root[a.max(b)] = a.min(b);
            }
        }
    }
    let mut order: Vec<(usize, u32)> = (0..n).map(|c| (find(&mut root, c), c as u32)).collect();
    order.sort_unstable();
    order.into_iter().map(|(_, c)| c).collect()
}

/// Clusters of up to `part` classes grown one class at a time: the next class is the unplaced
/// one with the most neighbours in the cluster (the smallest class among equals); a cluster
/// starts from the unplaced class with the most neighbours among the placed ones (the smallest
/// unplaced class when none has any). The clusters in the order grown, each in the order
/// grown.
pub fn grown_order(adj: &[Vec<u32>], part: usize) -> Vec<u32> {
    let n = adj.len();
    let part = part.max(1);
    let mut placed = vec![false; n];
    // neighbours in the current cluster
    let mut gain = vec![0u32; n];
    // neighbours among all placed classes (for the seeds)
    let mut seen = vec![0u32; n];
    let mut order: Vec<u32> = Vec::with_capacity(n);
    let best = |score: &[u32], placed: &[bool]| -> Option<usize> {
        (0..n)
            .filter(|&c| !placed[c])
            .max_by_key(|&c| (score[c], std::cmp::Reverse(c)))
    };
    while order.len() < n {
        let Some(seed) = best(&seen, &placed) else {
            break;
        };
        gain.fill(0);
        let mut size = 0;
        let mut next = Some(seed);
        while let Some(c) = next {
            placed[c] = true;
            order.push(c as u32);
            size += 1;
            for &d in &adj[c] {
                gain[d as usize] += 1;
                seen[d as usize] += 1;
            }
            if size >= part {
                break;
            }
            // the next class must touch the cluster: a cluster never spans components
            next = best(&gain, &placed).filter(|&c| gain[c] > 0);
        }
    }
    order
}

/// The class files a step reads to push into parts of `part` consecutive classes of `order`:
/// per part, the classes with a move into it (`adj` is symmetric: its neighbours), summed over
/// the parts.
pub fn reads(adj: &[Vec<u32>], order: &[u32], part: usize) -> usize {
    let n = adj.len();
    let mut mark = vec![usize::MAX; n];
    let mut total = 0;
    for (j, chunk) in order.chunks(part.max(1)).enumerate() {
        for &c in chunk {
            for &s in &adj[c as usize] {
                if mark[s as usize] != j {
                    mark[s as usize] = j;
                    total += 1;
                }
            }
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every order is a permutation of the classes, and the grown order reads fewer files
    /// than the class order on a cycle with chords.
    #[test]
    fn orders_are_permutations() {
        let n = 200usize;
        let mut adj: Vec<Vec<u32>> = vec![vec![]; n];
        for c in 0..n {
            for d in [(c + 1) % n, (c + 7) % n, (c + 50) % n] {
                adj[c].push(d as u32);
                adj[d].push(c as u32);
            }
        }
        for a in &mut adj {
            a.sort_unstable();
            a.dedup();
        }
        let grown = grown_order(&adj, 10);
        let mut sorted = grown.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..n as u32).collect::<Vec<_>>());
        let plain: Vec<u32> = (0..n as u32).collect();
        assert!(reads(&adj, &grown, 10) <= reads(&adj, &plain, 10));
        assert_eq!(reads(&adj, &plain, n), n);
    }
}
