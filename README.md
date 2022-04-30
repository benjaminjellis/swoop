[![CircleCI](https://circleci.com/gh/benjaminjellis/swoop/tree/master.svg?style=svg)](https://circleci.com/gh/benjaminjellis/swoop/tree/master)
![MSRV](https://img.shields.io/badge/msrv-1.60.0-red)
![version](https://img.shields.io/badge/version-0.1.0-yellow)
# swoop

Simple, lightweight optimisation algorithms in pure Rust 

## Motivation 

This crate aims to mimic the [scipy.optimize](https://docs.scipy.org/doc/scipy/reference/optimize.html) module in pure 
Rust.

## Example 
This crate has an asynchronous API and all examples use [Tokio](https://tokio.rs/).
To start your Cargo.toml should at least include


```toml
[dependencies]
swoop = { "git" = "https://github.com/benjaminjellis/swoop" }
tokio = { version = "1", features = ["full"] }
```

To minimise the function `f(x) = 3x^2 + 4x + 50` in the bound `-10 <= x <= 10` you can use the `bounded` optimiser

```rust
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
```