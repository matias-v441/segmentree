pub mod input{

use thiserror::Error;

#[derive(Debug, Error)]
pub enum InputError {
    #[error("Input contains NaN")]
    ContainsNaN,
    #[error("Input contains Infinite")]
    ContainsInfinite,
    #[error("Invalid interval: start > end")]
    StartGreaterThanEnd,
}
pub trait InputValidation{
    fn validate(self) -> Result<(), InputError>;
    fn validate_inf(self) -> Result<(), InputError>;
}

impl InputValidation for f64 {
    fn validate(self) -> Result<(), InputError> {
        if self.is_nan() {
            return Err(InputError::ContainsNaN);
        }
        if !self.is_finite() {
            return Err(InputError::ContainsInfinite);
        }
        Ok(())
    }
    fn validate_inf(self) -> Result<(), InputError> {
        if self.is_nan() {
            return Err(InputError::ContainsNaN);
        }
        Ok(())
    }
    
}

impl InputValidation for (f64, f64) {
    fn validate(self) -> Result<(), InputError> {
        self.validate_inf()?;
        if !self.0.is_finite() || !self.1.is_finite() {
            return Err(InputError::ContainsInfinite);
        }
        Ok(())
    }
    fn validate_inf(self) -> Result<(), InputError> {
        if self.0.is_nan() || self.1.is_nan() {
            return Err(InputError::ContainsNaN);
        }
        if self.1 < self.0 {
            return Err(InputError::StartGreaterThanEnd);
        }
        Ok(())
    }   
}

impl InputValidation for &[f64] {
    fn validate(self) -> Result<(), InputError> {
        for &x in self.iter() {
            x.validate()?;
        }
        Ok(())
    }
    fn validate_inf(self) -> Result<(), InputError> {
        for &x in self.iter() {
            x.validate_inf()?;
        }
        Ok(())
    }   
}

}