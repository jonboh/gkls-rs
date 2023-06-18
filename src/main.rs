use gkls_rs::{Problem, Options, ranf};
use plotly::{
    common::{ColorScale, ColorScalePalette, Marker, MarkerSymbol, Mode, Title},
    layout::{Axis, Layout},
    Mesh3D, Plot, Scatter3D, Surface, Configuration,
};




fn main() -> Result<(), Box<dyn std::error::Error>> {
    let problem = Problem::new(9, Options::default(), 2, 10, -1.0, 1./3., 2./3.).expect("Problem has to be valid");

    let n: i64 = 50;
    let xs: Vec<f64> = (-n..n).map(|f| f as f64 /n as f64).collect();
    let ys: Vec<f64> = (-n..n).map(|f| f as f64 /n as f64).collect();
    let zs: Vec<Vec<f64>> = xs
        .iter()
        .map(|x| {
            ys.iter()
                .map(|y| problem.d2_func(&[*x,*y]))
                .collect()
        })
        .collect();

    let trace = Surface::new(zs).x(xs).y(ys);
    let mut plot = Plot::new();
    plot.set_configuration(Configuration::new().autosizable(true).fill_frame(true));
    // plot.set_layout(Layout::new().auto_size(false).width(2000));
    plot.add_trace(trace);

    plot.show();
    Ok(())
}
#[test]
fn entry_point() {
    main().unwrap()
}
