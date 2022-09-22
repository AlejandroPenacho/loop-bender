use crate::app::diagrams::DiagramsConfiguration;
use crate::app::system::{Controller, DynamicalSystem, Model};

pub fn show_response_plot(
    ui: &mut egui::Ui,
    model: &Model,
    config: &DiagramsConfiguration,
    controller: &Controller,
    saved_controller: &Controller,
) {
    let system = model.link_system(controller);

    let (time, output) = compute_response(
        system,
        |_t, y| 1.0 - y,
        Timestep::Constant(0.002),
        ResponseDuration::Fixed(10.0),
    );

    let plot_line: egui::plot::PlotPoints = time
        .into_iter()
        .zip(output.into_iter())
        .map(|(t, y)| [t, y])
        .collect();

    let saved_controller_plot = if config.show_saved_controller {
        let system = model.link_system(saved_controller);
        let (time, output) = compute_response(
            system,
            |_t, y| 1.0 - y,
            Timestep::Constant(0.002),
            ResponseDuration::Fixed(10.0),
        );

        let plot_line: egui::plot::PlotPoints = time
            .into_iter()
            .zip(output.into_iter())
            .map(|(t, y)| [t, y])
            .collect();
        Some(plot_line)
    } else {
        None
    };

    egui::plot::Plot::new("response").show(ui, |ui| {
        ui.line(egui::plot::Line::new(plot_line));

        if let Some(x) = saved_controller_plot {
            ui.line(egui::plot::Line::new(x).style(egui::plot::LineStyle::Dashed { length: 1.5 }));
        }
    });
}

#[allow(dead_code)]
enum Timestep {
    Constant(f64),
    Automatic,
}

#[allow(dead_code)]
enum ResponseDuration {
    Fixed(f64),
    Automatic,
}

fn compute_response<T: DynamicalSystem>(
    system: T,
    input: fn(f64, f64) -> f64,
    timestep_conf: Timestep,
    duration: ResponseDuration,
) -> (Vec<f64>, Vec<f64>) {
    let mut time_vector = vec![0.0];
    let mut output_vector = vec![0.0];

    let mut model = system.get_state_space();
    let mut time = 0.0;

    let final_time = match duration {
        ResponseDuration::Fixed(x) => x,
        ResponseDuration::Automatic => unimplemented!(),
    };

    loop {
        let timestep = match timestep_conf {
            Timestep::Constant(x) => x,
            Timestep::Automatic => unimplemented!(),
        };

        let input = input(time, model.get_output());

        model.step(timestep, input);

        time += timestep;

        time_vector.push(time);
        output_vector.push(model.get_output());

        if time >= final_time {
            break;
        }
    }

    (time_vector, output_vector)
}
