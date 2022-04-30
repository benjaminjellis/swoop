use crate::minimise_scalar::{bracket_default, ScalarObjectiveFunction, ScalarOptimisationResult};
use crate::SwoopErrors;

/// Golden section univariate scalar optimisation
///
/// # Parameters
/// * `objective_function` - objective struct that implements the trait `ScalarObjectiveFunction`
/// * `xtol` - relative error in solution that is acceptable for convergence
/// * `maxiter` - maximum iterations
///
/// # Errors
/// Will return `SwoopErrors::MaxIterExceeded` if maximum number of iterations are exceeded before
/// the optimisation converges or `SwoopErrors::ArgumentError` is passed tolerance `xtol` is negative
pub async fn golden<T: ScalarObjectiveFunction>(
    objective_function: T,
    xtol: Option<f64>,
    maxiter: usize,
) -> Result<ScalarOptimisationResult, SwoopErrors> {
    let tol: f64;
    if let Some(i) = xtol {
        tol = i;
        if tol < 0f64 {
            return Err(SwoopErrors::ArgumentError(String::from(
                "Tolerance cannot be negative",
            )));
        }
    } else {
        tol = 2.22e-16;
    }

    let bracket = bracket_default(&objective_function)?;
    let (xa, xb, xc) = bracket.bracket;
    let mut fun_calls = bracket.fun_calls;

    // golden ratio conjugate: 2.0/(1.0+sqrt(5.0))
    let gr = 0.618_033_99;
    let gc = 1.0 - gr;
    let mut x3 = xc;
    let mut x0 = xa;
    let mut x1: f64;
    let mut x2: f64;

    if (xc - xb).abs() > (xb - xa).abs() {
        x1 = xb;
        x2 = xb + gc * (xc - xb);
    } else {
        x2 = xb;
        x1 = xb - gc * (xb - xa);
    }

    let mut f1 = objective_function.evaluate(x1);
    let mut f2 = objective_function.evaluate(x2);
    fun_calls += 1usize;

    let mut nit = 0usize;

    for _ in 0..maxiter {
        if (x3 - x0).abs() <= tol * (x1.abs() + x2.abs()) {
            break;
        }

        if f2 < f1 {
            x0 = x1;
            x1 = x2;
            x2 = gr * x1 + gc * x3;
            f1 = f2;
            f2 = objective_function.evaluate(x2);
        } else {
            x3 = x2;
            x2 = x1;
            x1 = gr * x2 + gc * x0;
            f2 = f1;
            f1 = objective_function.evaluate(x1);
        }

        fun_calls += 1;
        nit += 1;
    }

    let xmin: f64;
    let fval: f64;

    if f1 < f2 {
        xmin = x1;
        fval = f1;
    } else {
        xmin = x2;
        fval = f2;
    }

    let success: bool = nit < maxiter;

    if success {
        Ok(ScalarOptimisationResult {
            fun: fval,
            nfev: fun_calls,
            success: true,
            x: xmin,
        })
    } else {
        Err(SwoopErrors::MaxIterExceeded)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

    #[tokio::test]
    async fn test_quadratic() -> Result<(), SwoopErrors> {
        struct QuadraticFunction {}

        impl ScalarObjectiveFunction for QuadraticFunction {
            fn evaluate(&self, x: f64) -> f64 {
                (x - 2f64) * x * (x + 2f64).powf(2f64)
            }
        }
        let objective_function = QuadraticFunction {};
        let result = golden(objective_function, None, 500usize).await?;
        assert_eq!(result.success, true);
        assert_eq!(
            relative_eq!(result.fun, -9.914949590828147, epsilon = 1e-6),
            true
        );
        assert_eq!(
            relative_eq!(result.x, 1.2807764040333458, epsilon = 1e-6),
            true
        );
        Ok(())
    }
}
