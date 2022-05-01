use crate::minimise::{MultivariateObjectiveFunction, MultivariateOptimisationResult};
use crate::SwoopErrors;
use nalgebra::{clamp, DMatrix, DVector};


fn min_max(value: f64, min: f64, max: f64) -> f64{
    if value < min{
        min
    }else if value > max { max}
    else { value }
}

fn clip(mut initial_guess: DVector<f64>, lb: f64, ub: f64) -> DVector<f64> {
    initial_guess.map(|i| min_max(i, lb, ub))
}

/// Nelder Mead Simplex
///
/// # Parameters
/// todo
///
/// # References
///
/// [1] Gao, F. and Han, L.
///
///    Implementing the Nelder-Mead simplex algorithm with adaptive
///    parameters. 2012. Computational Optimization and Applications.
///    51:1, pp. 259-277
///
/// [2] [Scipy](https://github.com/scipy/scipy/blob/a6a2fe5e1f612aca080e2a150fd2a4c602ad10b6/scipy/optimize/_optimize.py#L635-L909)
///
///
pub async fn nelder_mead<T: MultivariateObjectiveFunction>(
    objective_functions: T,
    mut inital_guess: DVector<f64>,
    maxiter: usize,
    bounds: Option<(f64, f64)>,
    xtol: Option<f64>,
    fatol: Option<f64>,
) -> Result<MultivariateOptimisationResult, SwoopErrors> {
    let xtol_u: f64;
    let fatol_u: f64;

    if let Some(i) = xtol {
        xtol_u = i;
    } else {
        xtol_u = 1e-4f64;
    }

    if let Some(i) = xtol {
        fatol_u = i;
    } else {
        fatol_u = 1e-4f64;
    }

    let rho = 1f64;
    let chi = 2f64;
    let psi = 0.5f64;
    let sigma = 0.5f64;

    let nonzdelt = 0.05f64;
    let zdelt = 0.00025f64;


    if let Some(i) = bounds {
        if i.0 > i.1 {
            return Err(SwoopErrors::ArgumentError(String::from(
                "The lower bound exceeds the upper bound",
            )));
        }
        inital_guess = clip(inital_guess, i.0, i.1);
    }

    let n = inital_guess.len();
    let sim = DMatrix::


    Ok(MultivariateOptimisationResult {
        fun: 0.0,
        nfev: 0,
        success: false,
        x: DMatrix::from_vec(2, 1, vec![1.0, 1.0]),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_rosebrock() -> Result<(), SwoopErrors> {
        struct Rosenbrock {}
        impl MultivariateObjectiveFunction for Rosenbrock {
            fn evaluate(&self, x: DMatrix<f64>) -> f64 {
                let xx = x[(0, 0)];
                let y = x[(0, 0)];
                (1f64 - xx).powf(2f64) + 100f64 * (y - xx.powf(2f64)).powf(2f64)
            }
        }
        let rbf = Rosenbrock {};
        let ig = DVector::from_vec(vec![1.0, -20.0]);
        let bounds = (-10.0f64, 10.0f64);
        let resukt = nelder_mead(rbf, ig, 500, Some(bounds), None, None).await?;
        Ok(())
    }
}
