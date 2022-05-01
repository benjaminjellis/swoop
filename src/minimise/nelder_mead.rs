use crate::minimise::{MultivariateObjectiveFunction, MultivariateOptimisationResult};
use crate::SwoopErrors;
use nalgebra::{DMatrix};

struct Bounds {
    upper: DMatrix<f64>,
    lower: DMatrix<f64>,
}

/// Nelder Mead
///
/// References
///
/// [1] Gao, F. and Han, L.
///
///    Implementing the Nelder-Mead simplex algorithm with adaptive
///    parameters. 2012. Computational Optimization and Applications.
///    51:1, pp. 259-277
///
/// [2] [Scipy](https://github.com/scipy/scipy/blob/8fa4dcae4553fa873c5553e018f80a8a6f5e9948/scipy/optimize/_optimize.py#L635-L888)
///
///
pub async fn nelder_mead<T: MultivariateObjectiveFunction>(
    objective_functions: T,
    inital_guess: DMatrix<f64>,
    maxiter: usize,
    bounds: Option<DMatrix<f64>>,
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

    let formatted_bounds: Option<Bounds>;

    if let Some(i) = bounds {
        // todo why does column_sum() work here to pull out owned data?
        let m = i.column(0).column_sum().data;
        let t = m.as_vec();
        let n = i.column(1).column_sum().data;
        dbg!(t);
        dbg!(n);
    }

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
    async fn test_rosebrock() -> Result<(), SwoopErrors>{
        struct Rosenbrock {}
        impl MultivariateObjectiveFunction for Rosenbrock {
            fn evaluate(&self, x: DMatrix<f64>) -> f64 {
                let xx = x[(0, 0)];
                let y = x[(0, 0)];
                (1f64 - xx).powf(2f64) + 100f64 * (y - xx.powf(2f64)).powf(2f64)
            }
        }
        let rbf = Rosenbrock {};
        let ig = DMatrix::from_vec(2, 1, vec![1.0, 2.0]);
        let bounds =  DMatrix::from_vec(2, 2, vec![-10.0, -10.0, 10.0, 10.0]);
        let resukt = nelder_mead(rbf,
                                 ig,
            500,
            Some(bounds),
            None,
                                 None
        ).await?;
        Ok(())
    }
}
