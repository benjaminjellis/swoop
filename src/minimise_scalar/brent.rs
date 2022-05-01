use crate::minimise_scalar::{bracket_default, ScalarObjectiveFunction, ScalarOptimisationResult};
use crate::SwoopErrors;

/// Brent univariate scalar optimisation
///
/// # Parameters
/// * `objective_function` - objective struct that implements the trait `ScalarObjectiveFunction`
/// * `xtol` - relative error in solution that is acceptable for convergence
/// * `maxiter` - maximum iterations
///
/// # Errors
/// Will return `SwoopErrors::ArgumentError` is passed tolerance `xtol` is negative
///
/// # References
/// [1] [Scipy](https://github.com/scipy/scipy/blob/a6a2fe5e1f612aca080e2a150fd2a4c602ad10b6/scipy/optimize/_optimize.py#L2278-L2399)
#[allow(clippy::too_many_lines)]
pub async fn brent<T: ScalarObjectiveFunction>(
    objective_function: T,
    xtol: Option<f64>,
    maxiter: usize,
) -> Result<ScalarOptimisationResult, SwoopErrors> {
    let error_margin = f64::EPSILON;
    let tol: f64;
    if let Some(i) = xtol {
        tol = i;
        if tol < 0f64 {
            return Err(SwoopErrors::ArgumentError(String::from(
                "Tolerance cannot be negative",
            )));
        }
    } else {
        tol = 1.48e-8;
    }

    let bracket = bracket_default(&objective_function)?;
    let (xa, xb, xc) = bracket.bracket;
    //let (fa, fb, fc) = bracket.objective_function_values;
    let mut fun_calls = bracket.fun_calls;

    let min_tol = 1.0e-11;
    let cg = 0.381_966_0;

    let mut x = xb;
    let mut w = xb;
    let mut v = xb;
    let mut a: f64;
    let mut b: f64;

    let mut fw = objective_function.evaluate(x);
    let mut fv = fw;
    let mut fx = fw;

    if xa < xc {
        a = xa;
        b = xc;
    } else {
        a = xc;
        b = xa;
    }

    let mut deltax: f64 = 0.0;
    fun_calls += 1;
    let mut iter = 0usize;

    let mut tol1: f64;
    let mut tol2: f64;
    let mut xmid: f64;
    // is rat = 0.0f64 ok?
    let mut rat: f64 = deltax * cg;
    let mut tmp1: f64;
    let mut tmp2: f64;
    let mut p: f64;
    let mut dx_temp: f64;
    let mut u: f64;
    let mut fu: f64;

    while iter < maxiter {
        tol1 = tol * x.abs() + min_tol;
        tol2 = 2.0 * tol1;
        xmid = 0.5 * (a + b);

        // check for convergence
        if (x - xmid).abs() < (tol2 - 0.5 * (b - a)) {
            break;
        }

        if deltax.abs() <= tol1 {
            if x >= xmid {
                // do a golden section step
                deltax = a - x;
            } else {
                deltax = b - x;
            }
            rat = cg * deltax;
        } else {
            // do a parabolic step
            tmp1 = (x - w) * (fx - fv);
            tmp2 = (x - v) * (fx - fw);
            p = (x - v) * tmp2 - (x - w) * tmp1;
            tmp2 = 2.0 * (tmp2 - tmp1);
            if tmp2 > 0.0 {
                p = -p;
            }
            tmp2 = tmp2.abs();
            dx_temp = deltax;
            deltax = rat;
            // check parabolic fit
            if (p > tmp2 * (a - x))
                && (p < tmp2 * (b - x))
                && (p.abs() < (0.5 * tmp2 * dx_temp).abs())
            {
                rat = p * 1.0 / tmp2;
                u = x + rat;
                if (u - a) < tol2 || (b - u) < tol2 {
                    if xmid - x >= 0.0f64 {
                        rat = tol1;
                    } else {
                        rat = -tol1;
                    }
                }
            } else {
                if x >= xmid {
                    // if it's not do a golden section step
                    deltax = a - x;
                } else {
                    deltax = b - x;
                }
                rat = cg * deltax;
            }
        }
        if rat.abs() < tol1 {
            if rat >= 0.0f64 {
                u = x + tol1;
            } else {
                u = x - tol1;
            }
        } else {
            u = x + rat;
        }
        fu = objective_function.evaluate(u);
        fun_calls += 1;

        if fu > fx {
            if u < x {
                a = u;
            } else {
                b = u;
            }

            if (fu <= fw) || (w - x).abs() < error_margin {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if (fu <= fv) || (v - x).abs() < error_margin || (v - w).abs() < error_margin {
                v = u;
                fv = fu;
            }
        } else {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        iter += 1;
    }

    Ok(ScalarOptimisationResult {
        fun: fx,
        nfev: fun_calls,
        success: true,
        x,
    })
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
        let result = brent(objective_function, None, 500usize).await?;
        assert_eq!(result.success, true);
        assert_eq!(
            relative_eq!(result.fun, -9.914949590828147, epsilon = 1e-12),
            true
        );
        assert_eq!(
            relative_eq!(result.x, 1.2807764040333458, epsilon = 1e-12),
            true
        );
        Ok(())
    }
}
