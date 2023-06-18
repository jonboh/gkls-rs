use gkls_rs::ranf::Ranf;

fn main() {
    let mut rnd = Ranf::new(100, 37, 70, 1009, 1009, 200900);
    for _ in 0..2000 {
        println!("{:.6}", rnd.gen::<f64>());
    }
}
