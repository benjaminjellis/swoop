use swoop::minimise_scalar::{bounded, ScalarObjectiveFunction};
use swoop::SwoopErrors;

struct MyObjectiveFunction {
    a: f64,
    b: f64,
    c: f64,
}

impl MyObjectiveFunction {
    fn new(a: f64, b: f64, c: f64) -> Self {
        Self { a, b, c }
    }
}

impl ScalarObjectiveFunction for MyObjectiveFunction {
    fn evaluate(&self, x: f64) -> f64 {
        self.a * x.powf(2f64) + self.b * x + self.c
    }
}

#[tokio::main]
async fn main() -> Result<(), SwoopErrors> {
    let objective_function = MyObjectiveFunction::new(3f64, 4f64, 50f64);
    let result = bounded(objective_function, (-10f64, 10f64), 500usize).await?;
    println!("{:?}", result);
    Ok(())
}
