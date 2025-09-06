pub use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyDict;
use crate::core;

#[pyclass]
struct SegmentTree {
    inner: core::SegmentTree,
}

impl From<core::InputError> for PyErr {
    fn from(e: core::InputError) -> Self {
        PyValueError::new_err(e.to_string())
    }
}       

impl IntoPy<PyObject> for core::Stats {
    fn into_py(self, py: Python<'_>) -> PyObject {
        let d = PyDict::new_bound(py);
        d.set_item("length", self.length).unwrap();
        d.set_item("max_ovp", self.max_ovp).unwrap();
        d.set_item("min_ovp", self.min_ovp).unwrap();
        d.into_py(py)
    }
}

#[pyclass]
struct Union{
    inner: core::Union
}

#[pymethods]
impl Union {
    fn contains_point(&self, query: f64) -> bool {
        self.inner.contains_point(query)
    }
}

#[pymethods]
impl SegmentTree {

    #[new]
    fn new(mut all_ends: Vec<f64>) -> PyResult<Self> {
        Ok(Self {
            inner: core::SegmentTree::new(all_ends)?
        })
    }
    
    fn add_segment(&mut self, interval: (f64, f64), id: usize) -> PyResult<()> {
        self.inner.add_segment(interval, id)?;
        Ok(())
    }
    
    fn remove_segment(&mut self, interval: (f64, f64), id: usize) -> PyResult<()> {
        self.inner.remove_segment(interval, id)?;
        Ok(())
    }

    fn get_union(&self, interval: (f64, f64)) -> PyResult<Union> {
        Ok(Union { inner: self.inner.get_union(interval)? })
    }

    #[getter]
    fn root_stats(&self) -> core::Stats {
        self.inner.root_stats()
    }
}

pub fn register(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<SegmentTree>()?;
    m.add_class::<Union>()?;
    Ok(())
}
