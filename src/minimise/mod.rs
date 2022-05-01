//! Algorithms for multivariate function optimisation
use nalgebra as na;

mod nelder_mead;
mod slsqp;

/// Struct to represent the result of a scalar univariate function optimisation
#[derive(Debug, Clone)]
pub struct MultivariateOptimisationResult {
    /// Value of the objective function
    pub fun: f64,
    /// Number of evaluations of the objective function
    pub nfev: usize,
    /// Whether the optimisation was successful or not
    pub success: bool,
    /// The solution of the optimization
    pub x: na::DMatrix<f64>,
}

/// Trait to implement for a scalar multivariate objective function
pub trait MultivariateObjectiveFunction {
    /// Method to implement the objective function that will be used for evaluation when
    /// optimising
    fn evaluate(&self, x: na::DMatrix<f64>) -> f64;
}
