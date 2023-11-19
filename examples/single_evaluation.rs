use gkls_rs::{Options, Problem};
fn main() {
    let problem = Problem::new(9, Options::default(), 2, 10, -1.0, 1. / 3., 2. / 3.)
        .expect("Problem has to be valid");
    // Evaluate the proble with the appropiate function
    let x = [0.23, 0.44];
    let y = problem.d_func(&x);
    let dy = problem.d_gradient(&x);
    println!("f({x:?})={y}");
    println!("f'({x:?})={dy:?}");
}
