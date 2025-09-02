pub mod core;
mod util;

#[cfg(feature = "python")] pub mod py;
#[cfg(feature = "python")] use pyo3::prelude::*;
#[cfg(feature = "python")] 
#[pymodule]
fn segtree_native(_py: Python<'_>, m: &Bound<PyModule>) -> PyResult<()> {
   py::register(m)
}