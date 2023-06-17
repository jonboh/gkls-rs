use gkls_rs::{Problem, Options};
use plotters::prelude::*;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let problem = Problem::new(1, Options::default(), 2, 2, 0.0, 0.25, 0.25).unwrap();

    // PLOTTING
    let filepath = "plots/3d-plot.svg";
    let area = SVGBackend::new(filepath, (1024, 760)).into_drawing_area();

    area.fill(&WHITE)?;

    let x_axis = (-1.0..1.0).step(0.1);
    let z_axis = (-1.0..1.0).step(0.1);

    let mut chart = ChartBuilder::on(&area)
        .caption(format!("3D Plot Test"), ("sans", 20))
        .build_cartesian_3d(x_axis.clone(), -3.0..3.0, z_axis.clone())?;

    chart.with_projection(|mut pb| {
        pb.yaw = 0.5;
        pb.scale = 0.9;
        pb.into_matrix()
    });

    chart
        .configure_axes()
        .light_grid_style(BLACK.mix(0.15))
        .max_light_lines(3)
        .draw()?;

    chart
        .draw_series(
            SurfaceSeries::xoz(
                (-10..10).map(|f| f as f64 / 10.0),
                (-10..10).map(|f| f as f64 / 10.0),
                |x, y| (x * x + y * y).cos(),
            )
            .style(BLUE.mix(0.2).filled()),
        )?
        .label("Surface")
        .legend(|(x, y)| Rectangle::new([(x + 5, y - 5), (x + 15, y + 5)], BLUE.mix(0.5).filled()));
    // chart
    //     .draw_series(
    //         SurfaceSeries::xoz(
    //             (-10..10).map(|f| f as f64 / 10.0),
    //             (-10..10).map(|f| f as f64 / 10.0),
    //             |x, y| problem.nd_func(&[x,y]),
    //         )
    //         .style(BLUE.mix(0.2).filled()),
    //     )?
    //     .label("Surface")
    //     .legend(|(x, y)| Rectangle::new([(x + 5, y - 5), (x + 15, y + 5)], BLUE.mix(0.5).filled()));
    // chart
    //     .draw_series(
    //         SurfaceSeries::xoz(
    //             (-10..10).map(|f| f as f64 / 10.0),
    //             (-10..10).map(|f| f as f64 / 10.0),
    //             |x, y| problem.d_func(&[x,y]),
    //         )
    //         .style(BLUE.mix(0.2).filled()),
    //     )?
    //     .label("Surface")
    //     .legend(|(x, y)| Rectangle::new([(x + 5, y - 5), (x + 15, y + 5)], BLUE.mix(0.5).filled()));
    // chart
    //     .draw_series(
    //         SurfaceSeries::xoz(
    //             (-10..10).map(|f| f as f64 / 10.0),
    //             (-10..10).map(|f| f as f64 / 10.0),
    //             |x, y| problem.d2_func(&[x,y]),
    //         )
    //         .style(BLUE.mix(0.2).filled()),
    //     )?
    //     .label("Surface")
    //     .legend(|(x, y)| Rectangle::new([(x + 5, y - 5), (x + 15, y + 5)], BLUE.mix(0.5).filled()));
    //
    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .draw()?;

    // To avoid the IO failure being ignored silently, we manually call the present function
    area.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", filepath);
    Ok(())
}
#[test]
fn entry_point() {
    main().unwrap()
}
