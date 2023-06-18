#[link(name = "gkls", kind="static")]
extern "C" {
    fn nd_func(x: *const f64) -> f64;
}
