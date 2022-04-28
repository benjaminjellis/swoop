//! Algorithms for scalar univariate function optimisation

mod bounded;
pub use bounded::bounded;

mod brent;
pub use brent::brent;

mod golden;
use crate::SwoopErrors;
pub use golden::golden;

/// Struct to represent the result of a scalar univariate function optimisation
#[derive(Debug, Clone)]
pub struct ScalarOptimisationResult {
    /// Value of the objective function
    pub fun: f64,
    /// Number of evaluations of the objective function
    pub nfev: usize,
    /// Whether the optimisation was successful or not
    pub success: bool,
    /// The solution of the optimization
    pub x: f64,
}

/// Trait to implement for a scalar univarite objective function
pub trait ScalarObjectiveFunction {
    /// Method to implement the objective function that will be used for evaluation when
    /// optimising
    fn evaluate(&self, x: f64) -> f64;
}

pub(crate) struct BracketResult {
    pub(crate) bracket: (f64, f64, f64),
    pub(crate) fun_calls: usize,
}

/// Bracket the minimum of the function.
///
/// Given a function and distinct initial points, search in the
/// downhill direction (as defined by the initial points) and return
/// new points `xa, xb, xc` that bracket the minimum of the function
/// `f(xa) > f(xb) < f(xc)`.
///
/// It doesn't always mean that obtained solution will satisfy `xa<=x<=xb`.
///
/// # Parameters
/// * `objective_function` - Objective function to minimize
/// * `xa` - bottom of the bracketing interval
/// * `xb` - top of the bracketing interval
/// * `grow_limit` - Maximum grow limit
/// * `maxiter` - Maximum number of iterations to perform.
pub(crate) fn bracket<T: ScalarObjectiveFunction>(
    objective_function: &T,
    mut xa: f64,
    mut xb: f64,
    grow_limit: f64,
    maxiter: usize,
) -> Result<BracketResult, SwoopErrors> {
    // the golden ratio
    let gold = 1.618_034;

    let very_small_number = 1e-21;
    let mut fa = objective_function.evaluate(xa);
    let mut fb = objective_function.evaluate(xb);

    if fa < fb {
        (xa, xb) = (xb, xa);
        (fa, fb) = (fb, fa);
    }

    let mut xc = xb + gold * (xb - xa);
    let mut fc = objective_function.evaluate(xc);
    let mut fun_calls = 3usize;
    let mut iter = 0usize;

    let mut tmp1: f64;
    let mut tmp2: f64;
    let mut val: f64;

    let mut denom: f64;
    let mut w: f64;
    let mut wlim: f64;

    let mut fw: f64;

    while fc < fb {
        tmp1 = (xb - xa) * (fb - fc);
        tmp2 = (xb - xc) * (fb - fa);
        val = tmp2 - tmp1;

        if val.abs() < very_small_number {
            denom = 2.0 * very_small_number;
        } else {
            denom = 2.0 * val;
        }

        w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom;
        wlim = xb + grow_limit * (xc - xb);

        if iter > maxiter {
            return Err(SwoopErrors::MaxIterExceeded);
        }
        iter += 1usize;

        if (w - xc) * (xb - w) > 0.0 {
            fw = objective_function.evaluate(w);
            fun_calls += 1usize;
            if fw < fc {
                xa = xb;
                xb = w;
                //fa = fb;
                //fb = fw;
                return Ok(BracketResult {
                    bracket: (xa, xb, xc),
                    fun_calls,
                });
            } else if fw > fb {
                xc = w;
                //fc = fw;
                return Ok(BracketResult {
                    bracket: (xa, xb, xc),
                    fun_calls,
                });
            }

            w = xc + gold * (xc - xb);
            fw = objective_function.evaluate(w);
            fun_calls += 1usize;
        } else if (w - wlim) * (wlim - xc) >= 0.0 {
            w = wlim;
            fw = objective_function.evaluate(w);
            fun_calls += 1usize;
        } else if (w - wlim) * (xc - w) > 0.0 {
            fw = objective_function.evaluate(w);
            fun_calls += 1usize;
            if fw < fc {
                xb = xc;
                xc = w;
                w = xc + gold * (xc - xb);
                fb = fc;
                fc = fw;
                fw = objective_function.evaluate(w);
                fun_calls += 1usize;
            }
        } else {
            w = xc + gold * (xc - xb);
            fw = objective_function.evaluate(w);
            fun_calls += 1usize;
        }
        xa = xb;
        xb = xc;
        xc = w;
        fa = fb;
        fb = fc;
        fc = fw;
    }

    Ok(BracketResult {
        bracket: (xa, xb, xc),
        fun_calls,
    })
}

/// Bracket the minimum of the function with default arguments
///
/// # Parameters
/// * `objective_function` - Objective function to minimize
pub(crate) fn bracket_default<T: ScalarObjectiveFunction>(
    objective_function: &T,
) -> Result<BracketResult, SwoopErrors> {
    bracket(objective_function, 0.0f64, 1.0f64, 110.0f64, 1000usize)
}
