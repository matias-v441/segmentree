use core::f64;
use std::vec;
use thiserror::Error;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyDict;

#[derive(Default, Clone)]
struct Node {
    segs: Vec<usize>,
    itv: (f64, f64)
}

#[pyclass]
pub struct SegmentTree {
    nodes: Vec<Node>,
    stats: Vec<Stats>
}

trait ChildIdUnchecked { fn left(self) -> Self; fn right(self) -> Self; }
impl ChildIdUnchecked for usize {
    #[inline] fn left(self) -> Self { (self << 1) | 1 }
    #[inline] fn right(self) -> Self { (self << 1) + 2 }
}

#[derive(Debug, Error)]
pub enum IntervalError {
    #[error("Interval contains NaN")]
    ContainsNaN,
    #[error("Interval contains Infinite")]
    ContainsInfinite,
    #[error("Invalid interval: start > end")]
    StartGreaterThanEnd,
}

impl From<IntervalError> for PyErr {
    fn from(e: IntervalError) -> Self {
        PyValueError::new_err(e.to_string())
    }
}       

trait Interval {
    fn overlaps(self, other: Self) -> bool;
    fn contains(self, other: Self) -> bool;
    fn validate(self) -> Result<(), IntervalError>;
    fn validate_inf(self) -> Result<(), IntervalError>;
}

impl Interval for (f64, f64) {
    #[inline]
    fn overlaps(self, other: Self) -> bool {
        self.0 < other.1 && other.0 < self.1
    }
    fn contains(self, other: Self) -> bool {
        self.0 <= other.0 && other.1 <= self.1
    }
    fn validate(self) -> Result<(), IntervalError> {
        self.validate_inf()?;
        if !self.0.is_finite() || !self.1.is_finite() {
            return Err(IntervalError::ContainsInfinite);
        }
        Ok(())
    }
    fn validate_inf(self) -> Result<(), IntervalError> {
        if self.0.is_nan() || self.1.is_nan() {
            return Err(IntervalError::ContainsNaN);
        }
        if self.1 < self.0 {
            return Err(IntervalError::StartGreaterThanEnd);
        }
        Ok(())
    }
}


#[derive(Default, Clone, Copy)]
struct Stats {
    length: f64,
    max_ovp: usize,
}

impl IntoPy<PyObject> for Stats {
    fn into_py(self, py: Python<'_>) -> PyObject {
        let d = PyDict::new_bound(py);
        d.set_item("length", self.length).unwrap();
        d.set_item("max_ovp", self.max_ovp).unwrap();
        d.into_py(py)
    }
}

#[pyclass]
struct Union{
    intervals: Vec<(f64, f64)>
}

#[pymethods]
impl Union {
    fn contains_point(&self, query: f64) -> bool {
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

#[pymethods]
impl SegmentTree {

    #[new]
    fn new(mut all_ends: Vec<f64>) -> PyResult<Self> {
        if let Some(i) = all_ends.iter().position(|&x| x.is_nan() || x.is_infinite()) {
            return Err(PyErr::new::<PyValueError, _>(
                format!("Invalid value at index {}: NaN or Infinite", i),
            ));
        }
        let mut tree = Self { nodes: Vec::new(), stats: Vec::new() };
        all_ends.sort_by(|a, b| a.partial_cmp(b).unwrap());
        all_ends.dedup();
        tree.build(&all_ends);
        tree.stats = vec![Stats::default(); tree.nodes.len()];
        Ok(tree)
    }
    
    fn add_segment(&mut self, interval: (f64, f64), id: usize) -> PyResult<()> {
        interval.validate()?;
        self.seg_nodes_apply(interval,
            &mut |s, i| {
                s.nodes[i].segs.push(id)
            }, 0);
        Ok(())
    }
    
    fn remove_segment(&mut self, interval: (f64, f64), id: usize) -> PyResult<()> {
        interval.validate()?;
        self.seg_nodes_apply(interval, &mut |s, i| {
            s.nodes[i].segs.retain(|&x| x != id)
        }, 0);
        Ok(())
    }

    fn get_union(&self, interval: (f64, f64)) -> PyResult<Union> {
        interval.validate()?;
        let mut union = Union { intervals: Vec::new() };
        self.seg_nodes_visit(interval, &mut |s, i| {
            let node = &s.nodes[i];
            union.intervals.push(node.itv);
        }, 0);
        Ok(union)
    }

    #[getter]
    fn root_stats(&self) -> Stats {
        self.stats[0]
    }
}


#[pymodule]
fn segtree_native(_py: Python<'_>, m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<SegmentTree>()?;
    m.add_class::<Union>()?;
    Ok(())
}