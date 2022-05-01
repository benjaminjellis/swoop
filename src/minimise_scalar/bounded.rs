use crate::minimise_scalar::{ScalarObjectiveFunction, ScalarOptimisationResult};
use crate::SwoopErrors;

/// Returns -1 if x < 0, 0 if x==0, 1 if x > 0. nan is returned for nan inputs.
fn sign(number: f64) -> f64 {
    if (number - 0f64).abs() < f64::EPSILON {
        0f64
    } else if number > 0f64 {
        1f64
    } else {
        -1f64
    }
}

/// Check is a - b is zero, if it is return 1, else 0
fn zero_diff(a: f64, b: f64) -> f64 {
    if ((a - b) - 0f64).abs() < f64::EPSILON {
        1f64
    } else {
        0f64
    }
}

/// Check if a value is 0, if it is return 1 else return 0
fn zero_or_not(a: f64) -> f64 {
    if (a - 0f64).abs() < f64::EPSILON {
        1f64
    } else {
        0f64
    }
}

/// Find the maximum of the two arguments
///
/// # Parameters
/// * `a` - arg 1
/// * `b` - arg 1
fn arg_max(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else if b > a {
        b
    } else {
        a
    }
}

/// Bounded univariate scalar optimisation
///
/// # Parameters
/// * `objective_function` - objective struct that implements the trait `ScalarObjectiveFunction`
/// * `bounds` - bounds for the optimisation
/// * `maxiter` - maximum iterations
///
/// # Errors
/// Will return `SwoopErrors::ArgumentError` is passed tolerance `xtol` is negative
///
/// # References
/// [1] [Scipy](https://github.com/scipy/scipy/blob/a6a2fe5e1f612aca080e2a150fd2a4c602ad10b6/scipy/optimize/_optimize.py#L2082-L2224)
///
#[allow(clippy::cast_sign_loss)]
#[allow(clippy::cast_possible_truncation)]
#[allow(clippy::cast_precision_loss)]
#[allow(clippy::too_many_lines)]
pub async fn bounded<T: ScalarObjectiveFunction>(
    objective_function: T,
    bounds: (f64, f64),
    maxiter: usize,
) -> Result<ScalarOptimisationResult, SwoopErrors> {
    let error_margin = f64::EPSILON;
    let xatol = 1e-5f64;
    let (x1, x2) = bounds;
    if x2 < x1 {
        return Err(SwoopErrors::ArgumentError(String::from(
            "The lower bound exceeds the upper bound",
        )));
    }
    let sqrt_eps = f64::EPSILON.sqrt();
    let golden_mean = 0.5f64 * (3.0f64 - (5.0f64.sqrt()));
    let (mut a, mut b) = (x1, x2);
    let mut fulc = a + golden_mean * (b - a);
    let (mut nfc, mut xf) = (fulc, fulc);

    let mut rat = 0.0f64;
    let mut e = 0.0f64;

    let mut x = xf;

    let mut fx = objective_function.evaluate(x);
    let mut num = 1f64;

    let mut fu: f64;

    let mut ffulc = fx;
    let mut fnfc = fx;

    let mut xm = 0.5f64 * (a + b);

    let mut tol1 = sqrt_eps * (xf).abs() + xatol / 3.0f64;
    let mut tol2 = 2.0f64 * tol1;

    let mut golden: bool;
    let mut r: f64;
    let mut q: f64;
    let mut p: f64;
    let mut si: f64;

    let mut success = true;

    while (xf - xm).abs() > (tol2 - 0.5 * (b - a)) {
        golden = true;

        // Check for parabolic fit
        if e.abs() > tol1 {
            golden = true;
            r = (xf - nfc) * (fx - ffulc);
            q = (xf - ffulc) * (fx - fnfc);
            p = (xf - fulc) * q - (xf - nfc) * r;
            q = 2.0 * (q - r);
            if q > 0.0f64 {
                p = -p;
            }
            q = q.abs();
            r = e;
            e = rat;

            // Check for acceptability of parabola
            if (p.abs() < (0.5f64 * q * r).abs()) && (p > q * (a - xf)) && (p < q * (b - xf)) {
                rat = (p + 0.0f64) / q;
                x = xf + rat;

                if ((x - a) < tol2) || ((b - x) < tol2) {
                    si = sign(xm - xf) + zero_diff(xm, xf);
                    rat = tol1 * si;
                }
            } else {
                golden = true;
            }
        }

        if golden {
            if xf >= xm {
                e = a - xf;
            } else {
                e = b - xf;
            }

            rat = golden_mean * e;
        }

        si = sign(rat) + zero_or_not(rat);
        x = xf + si * arg_max(rat.abs(), tol1);
        fu = objective_function.evaluate(x);
        num += 1f64;

        if fu <= fx {
            if x >= xf {
                a = xf;
            } else {
                b = xf;
            }

            (fulc, ffulc) = (nfc, fnfc);
            (nfc, fnfc) = (xf, fx);
            (xf, fx) = (x, fu);
        } else {
            if x < xf {
                a = x;
            } else {
                b = x;
            }

            if (fu <= fnfc) || (nfc - xf).abs() < error_margin {
                (fulc, ffulc) = (nfc, fnfc);
                (nfc, fnfc) = (x, fu);
            } else if (fu <= ffulc)
                || (fulc - xf).abs() < error_margin
                || (fulc - nfc).abs() < error_margin
            {
                (fulc, ffulc) = (x, fu);
            }
        }

        xm = 0.5f64 * (a + b);
        tol1 = sqrt_eps * xf.abs() + xatol / 3.0f64;
        tol2 = 2.0 * tol1;

        if num >= maxiter as f64 {
            success = false;
            break;
        }
    }

    Ok(ScalarOptimisationResult {
        fun: fx,
        nfev: num as usize,
        success,
        x: xf,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::relative_eq;

    #[tokio::test]
    async fn test_quadratic() -> Result<(), SwoopErrors> {
        struct QuadraticFunction {
            a: f64,
            b: f64,
            c: f64,
        }

        impl QuadraticFunction {
            fn new(a: f64, b: f64, c: f64) -> Self {
                Self { a, b, c }
            }
        }

        impl ScalarObjectiveFunction for QuadraticFunction {
            fn evaluate(&self, x: f64) -> f64 {
                self.a * x.powf(2f64) + self.b * x + self.c
            }
        }

        let objective_function = QuadraticFunction::new(3f64, 4f64, 50f64);
        let result = bounded(objective_function, (-10f64, 10f64), 500usize).await?;
        println!("{:?}", result);
        assert_eq!(
            relative_eq!(result.fun, 48.666666666666664, epsilon = 1e-6),
            true
        );
        assert_eq!(
            relative_eq!(result.x, -0.666666666666667, epsilon = 1e-6),
            true
        );
        Ok(())
    }
}
