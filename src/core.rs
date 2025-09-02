use core::f64;
use std::vec;
pub use crate::util::input::*;

#[derive(Default, Clone)]
pub struct Node {
    segs: Vec<usize>,
    itv: (f64, f64)
}

pub struct SegmentTree {
    nodes: Vec<Node>,
    stats: Vec<Stats>
}

trait ChildIdUnchecked { fn left(self) -> Self; fn right(self) -> Self; }
impl ChildIdUnchecked for usize {
    #[inline] fn left(self) -> Self { (self << 1) | 1 }
    #[inline] fn right(self) -> Self { (self << 1) + 2 }
}


trait Interval {
    fn overlaps(self, other: Self) -> bool;
    fn contains(self, other: Self) -> bool;
}

impl Interval for (f64, f64) {
    #[inline]
    fn overlaps(self, other: Self) -> bool {
        self.0 < other.1 && other.0 < self.1
    }
    fn contains(self, other: Self) -> bool {
        self.0 <= other.0 && other.1 <= self.1
    }
}

#[derive(Default, Clone, Copy)]
pub struct Stats {
    pub length: f64,
    pub max_ovp: usize,
}

pub struct Union{
    pub intervals: Vec<(f64, f64)>
}

impl Union {
    pub fn contains_point(&self, query: f64) -> bool {
        if query.is_nan() || query.is_infinite() {
            return false;
        }
        let pos = self.intervals.binary_search_by(|itv| {
            if itv.0 > query {
                std::cmp::Ordering::Greater
            } else if itv.1 < query {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Equal
            }
        });
        pos.is_ok()
    }
}

impl SegmentTree {

    fn build(&mut self, ends: &[f64]) {
        let n_leaves = ends.len() << 1 + 1;
        let n_nodes_compl = n_leaves.next_power_of_two() >> 1; // nodes in last complete level
        let n_leaves_last = (n_leaves - n_nodes_compl) << 1; // leaves in last level
        let n_leaves_compl = n_leaves - n_leaves_last; // leaves in last complete level
        let n_nodes = n_nodes_compl << 1 - 1 + n_leaves_last;
        self.nodes = vec![Node::default(); n_nodes];
        let p_last_start = n_nodes_compl << 1 - 1;
        let p_leave_compl_start = p_last_start - n_leaves_compl - 1;
        self.nodes[p_last_start].itv.0 = f64::NEG_INFINITY;
        let leaf_ids:Vec<usize> = (p_last_start..n_nodes).chain(p_leave_compl_start..p_last_start).collect();
        for (i,&v) in ends.iter().enumerate() {
            let seg_prev = &mut self.nodes[leaf_ids[i << 1]];
            seg_prev.itv.1 = v;
            let point = &mut self.nodes[leaf_ids[i << 1 + 1]];
            point.itv.0 = v;
            point.itv.1 = v;
            let seg_next = &mut self.nodes[leaf_ids[i << 1 + 2]];
            seg_next.itv.0 = v;
        }
        self.nodes[n_nodes-1].itv.1 = f64::INFINITY;
        for i in (n_nodes-n_leaves-1)..0 {
            self.nodes[i].itv.0 = self.nodes[i.left()].itv.0;
            self.nodes[i].itv.1 = self.nodes[i.right()].itv.1;
        }
    }

    fn seg_nodes_apply(&mut self, seg: (f64, f64), op: &mut dyn FnMut(&mut Self, usize), node_id: usize) {
        if seg.contains(self.nodes[node_id].itv) {
            op(self, node_id);
            self.update_stats(node_id);
            return;
        }
        if self.is_leaf(node_id) {
            return;
        }
        if self.nodes[node_id.left()].itv.overlaps(seg){
            self.seg_nodes_apply(seg, op, node_id.left());
        }
        if self.nodes[node_id.right()].itv.overlaps(seg){
            self.seg_nodes_apply(seg, op, node_id.right());
        }
        self.update_stats(node_id);
    }

    fn seg_nodes_visit(&self, seg: (f64, f64), report: &mut dyn FnMut(&Self, usize), node_id: usize) {
        if seg.contains(self.nodes[node_id].itv) {
            report(self, node_id);
            return;
        }
        if self.is_leaf(node_id) {
            return;
        }
        if self.nodes[node_id.left()].itv.overlaps(seg){
            self.seg_nodes_visit(seg, report, node_id.left());
        }
        if self.nodes[node_id.right()].itv.overlaps(seg){
            self.seg_nodes_visit(seg, report, node_id.right());
        }
    }

    fn update_stats(&mut self, node_id: usize) {
        let node: &Node = &self.nodes[node_id];
        let mut new_stats = Stats {
            max_ovp: node.segs.len(),
            length: node.itv.1 - node.itv.0
        };
        if !self.is_leaf(node_id) {
            let child_stats = (self.stats[node_id.left()], self.stats[node_id.right()]);
            new_stats.max_ovp += Ord::max(child_stats.0.max_ovp, child_stats.1.max_ovp);
            if node.segs.is_empty() {
                new_stats.length = child_stats.0.length + child_stats.1.length;
            }
        }
    }

    #[inline]
    fn is_leaf(&self, idx: usize) -> bool {
        idx >= (self.nodes.len() >> 1)
    }
}

impl SegmentTree {

    pub fn new(mut all_ends: Vec<f64>) -> Result<Self, InputError> {
        all_ends.validate()?;
        let mut tree = Self { nodes: Vec::new(), stats: Vec::new() };
        all_ends.sort_by(|a, b| a.partial_cmp(b).unwrap());
        all_ends.dedup();
        tree.build(&all_ends);
        tree.stats = vec![Stats::default(); tree.nodes.len()];
        Ok(tree)
    }

    pub fn add_segment(&mut self, interval: (f64, f64), id: usize) -> Result<(), InputError> {
        interval.validate()?;
        self.seg_nodes_apply(interval,
            &mut |s, i| {
                s.nodes[i].segs.push(id)
            }, 0);
        Ok(())
    }

    pub fn remove_segment(&mut self, interval: (f64, f64), id: usize) -> Result<(), InputError> {
        interval.validate()?;
        self.seg_nodes_apply(interval, &mut |s, i| {
            s.nodes[i].segs.retain(|&x| x != id)
        }, 0);
        Ok(())
    }

    pub fn get_union(&self, interval: (f64, f64)) -> Result<Union, InputError> {
        interval.validate_inf()?;
        let mut union = Union { intervals: Vec::new() };
        self.seg_nodes_visit(interval, &mut |s, i| {
            let node = &s.nodes[i];
            union.intervals.push(node.itv);
        }, 0);
        Ok(union)
    }

    pub fn root_stats(&self) -> Stats {
        self.stats[0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_segment_tree() {
        let mut tree = SegmentTree::new(vec![1.0, 2.0, 3.0]).unwrap();
        tree.add_segment((1.0, 2.0), 0).unwrap();
        tree.add_segment((2.0, 3.0), 1).unwrap();
        assert_eq!(tree.root_stats().max_ovp, 2);
    }
}
