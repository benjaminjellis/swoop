//! Simple optimisation algorithms in pure rust
#![warn(clippy::pedantic)]
#![warn(missing_docs)]
#![allow(clippy::similar_names)]
#![allow(clippy::unused_async)]
#![allow(clippy::many_single_char_names)]

pub mod minimise;
pub mod minimise_scalar;

use thiserror::Error;

/// Error enum
#[derive(Error, Debug)]
pub enum SwoopErrors {
    /// Error to catch maximum iteration exceptions
    #[error("Maximum number of iterations exceeded")]
    MaxIterExceeded,
    /// Error for incorrectly set argument
    #[error("Invalid argument received `{0}`")]
    ArgumentError(String),
    /// Transparent error handler
    #[error(transparent)]
    Other(#[from] anyhow::Error), // source and Display delegate to anyhow::Error
}
