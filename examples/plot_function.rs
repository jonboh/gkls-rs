use plotters::prelude::*;

use gkls_rs::{Options, Problem};

const OUT_FOLDER: &str = "example_output";
const NF: usize = 9;

fn min_max(numbers: &[f64]) -> Option<(f64, f64)> {
    numbers.iter().fold(None, |min_max, &val| match min_max {
        None => Some((val, val)),
        Some((min, max)) => Some((min.min(val), max.max(val))),
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    std::fs::create_dir_all(&OUT_FOLDER).expect("Could not create output folder");
    let file_name = format!("{OUT_FOLDER}/nf{NF}.svg");
    let area = SVGBackend::new(&file_name, (1024, 760)).into_drawing_area();
    // let file_name = format!("{OUT_FOLDER}/nf{NF}.png");
    // let area = BitMapBackend::new(&file_name, (1024, 760)).into_drawing_area();

    let problem = Problem::new(9, Options::default(), 2, 10, -1.0, 1. / 3., 2. / 3.)
        .expect("Problem has to be valid");
    let n: i64 = 50;
    let xs: Vec<f64> = (-n..n).map(|f| f as f64 / n as f64).collect();
    let ys: Vec<f64> = (-n..n).map(|f| f as f64 / n as f64).collect();
    let zs: Vec<Vec<f64>> = xs
        .iter()
        .map(|x| ys.iter().map(|y| problem.d2_func(&[*x, *y])).collect())
        .collect();
    let zflatten: Vec<f64> = zs.iter().flat_map(|num| num.clone()).collect();
    let (xmin, xmax) = min_max(&xs).expect("X values should not have NaN");
    let (ymin, ymax) = min_max(&ys).expect("Y values should not have NaN");
    let (zmin, zmax) = min_max(&zflatten).expect("Z values should not have NaN");

    area.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&area)
        .caption(format!("Function N{NF}"), ("sans", 20))
        .build_cartesian_3d(xmin..xmax, zmin..zmax, ymax..ymin)?;

    chart.set_3d_pixel_range((500, 300, 500));

    chart.with_projection(|mut pb| {
        pb.yaw = -1.0;
        pb.pitch = 0.35;
        pb.scale = 1.15;
        pb.into_matrix()
    });

    chart.configure_axes().max_light_lines(3).draw()?;

    chart
        .draw_series(
            SurfaceSeries::xoz(xs.into_iter(), ys.into_iter(), |x, y| {
                problem.d_func(&[x, y])
            })
            .style_func(&|&v| {
                (ViridisRGB::get_color_normalized(v * 1.5, zmin, zmax))
                    .filled()
                    .into()
            }),
        )?
        .label(format!("Function N{NF} D-type"));

    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {file_name}");
    Ok(())
}
