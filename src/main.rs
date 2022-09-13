mod complex;
mod system;

use system::{PID, Model, DynamicalSystem};

use egui::widgets::Widget;

fn main() {
    eframe::run_native("My first App", eframe::NativeOptions::default(), Box::new(|cc| Box::new(MyApp::new(cc))));
}

struct MyApp {
    model: Model,
    pointer_mode: PointerMode,
    pid: PID,
    bode_link: egui::widgets::plot::LinkedAxisGroup
}

#[derive(PartialEq,Eq)]
enum PointerMode {
    AddPole,
    AddZero,
    Remove,
    Move
}

impl PointerMode {
    fn display(&self) -> &'static str {
       match self {
            PointerMode::AddPole => "Pole",
            PointerMode::AddZero => "Zero",
            PointerMode::Remove => "Remove",
            PointerMode::Move => "Move"
        }
    }
}

impl MyApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            model: Model::default(),
            pointer_mode: PointerMode::AddPole,
            pid: PID::default(),
            bode_link: egui::widgets::plot::LinkedAxisGroup::new(true, false)
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::containers::Window::new("Diagrams")
            .default_pos((0.0,0.0))
            .show(ctx, |ui| {

            ui.heading("Heyyy");
            get_bode_plot(ui, &self.model, &self.pid, &self.bode_link);
        });

        egui::containers::Window::new("Tuner")
            .default_pos((100.0,0.0))
            .show(ctx, |ui| {

            get_pz_map(ui, &mut self.pointer_mode, &mut self.model, &mut self.pid);
            // ui.add(egui::widgets::Label::new(format!("Poles: {:?}", self.model.get_poles())));
            // ui.add(egui::widgets::Label::new(format!("Zeros: {:?}", self.model.get_zeros())));
        });

        egui::containers::Window::new("Response")
            .default_pos((100.0,30.0))
            .show(ctx, |ui| {
                response_plot(ui, &self.model, &self.pid);
            });
    }
}

fn response_plot(ui: &mut egui::Ui, model: &Model, pid: &PID) {
    let mut state_space = model.get_state_space();

    let mut output: Vec<f64> = Vec::with_capacity(5000);
    let time = (0..5000).map(|x| x as f64/100.0).collect::<Vec<f64>>();

    output.push(0.0);
    let k = pid.borrow_k();

    for _ in 0..1000 {
        output.push(state_space.step(0.01, 1.0 - k * output.iter().last().unwrap()));//, 1.0- k * output.iter().last().unwrap()));
        // output.push(state_space.step(0.01, 1.0));//, 1.0- k * output.iter().last().unwrap()));
    }

    let plot_line: egui::plot::PlotPoints = time.into_iter().zip(output.into_iter()).map(|(t,y)| [t,y]).collect();

    egui::plot::Plot::new("response")
        .show(ui, |ui| {
            ui.line(egui::plot::Line::new(plot_line));
        });
}

fn pid_config(ui: &mut egui::Ui, pid: &mut PID) {
    ui.horizontal(|ui| {
        ui.vertical(|ui| {
            ui.horizontal(|ui| {
                ui.add_visible(false, egui::widgets::Checkbox::new(&mut false, ""));
                ui.add(egui::widgets::Label::new("K"));
            });
            ui.horizontal(|ui| {
                ui.add(egui::widgets::Checkbox::new(pid.get_integral_on(), ""));
                ui.add(egui::widgets::Label::new("T_I"));
            });
            ui.horizontal(|ui| {
                ui.add(egui::widgets::Checkbox::new(pid.get_derivative_on(), ""));
                ui.add(egui::widgets::Label::new("T_D"));
            });
        });

        let k_speed = *pid.get_k()/20.0;
        let i_speed = *pid.get_ti()/20.0;
        let d_speed = *pid.get_td()/20.0;

        ui.vertical(|ui| {
            ui.add(
                egui::widgets::DragValue::new(pid.get_k())
                    .clamp_range(0.01..=50.0)
                    .speed(k_speed)
            );

            ui.add_enabled(
                *pid.get_integral_on(),
                egui::widgets::DragValue::new(pid.get_ti())
                    .clamp_range(0.1..=50.0)
                    .speed(i_speed)
            );

            ui.add_enabled(
                *pid.get_derivative_on(),
                egui::widgets::DragValue::new(pid.get_td())
                    .clamp_range(0.1..=50.0)
                    .speed(d_speed)
            );
        });
    });
}

fn get_bode_plot(ui: &mut egui::Ui, model: &Model, pid: &PID, bode_link: &egui::widgets::plot::LinkedAxisGroup) {
    let all_freqs_expo = (0..500).map(|i| (i as f64/100.0) - 3.0).collect::<Vec<f64>>();
    let all_freqs = all_freqs_expo.iter().map(|i| 10f64.powf(*i)).collect::<Vec<f64>>();

    let mags = all_freqs.iter().map(|f| {
        let mut mag = 1f64;

        mag *= model.get_mag_at_freq(*f);
        mag *= pid.get_mag_at_freq(*f);

        mag.log10()

    }).collect::<Vec<f64>>();
    
    let phases = all_freqs.iter().map(|f| {
        let mut phase = 0.0;

        phase += model.get_phase_at_freq(*f);
        phase += pid.get_phase_at_freq(*f);

        phase

    }).collect::<Vec<f64>>();

    let mag_line: egui::plot::PlotPoints = (0..500).map(|i| [all_freqs_expo[i], mags[i]]).collect();
    let phase_line: egui::plot::PlotPoints  = (0..500).map(|i| [all_freqs_expo[i], phases[i]]).collect();

    let nyquist_line: egui::plot::PlotPoints =
        (0..500).map(|i| [10f64.powf(mags[i])*phases[i].cos(), 10f64.powf(mags[i])*phases[i].sin()])
        .collect();


    ui.vertical(|ui| {
        // ui.add(egui::widgets::Label::new(format!("{:?}", mags)));
        egui::plot::Plot::new("bode_mag")
            .view_aspect(2.0)
            .center_y_axis(false)
            .link_axis(bode_link.clone())
            .show(ui,
            |ui| {
                ui.line(egui::plot::Line::new(mag_line));
            });
        egui::plot::Plot::new("bode_phase")
            .view_aspect(2.0)
            .center_y_axis(false)
            .link_axis(bode_link.clone())
            .show(ui,
            |ui| {
                ui.line(egui::plot::Line::new(phase_line));
            });
        egui::plot::Plot::new("nyquist")
            .view_aspect(1.0)
            .data_aspect(1.0)
            .include_x(0.0)
            .include_y(0.0)
            .show(ui,
            |ui| {
                ui.line(egui::plot::Line::new(nyquist_line));
            });
    });
}

fn get_pz_map(
        ui: &mut egui::Ui,
        pointer_mode: &mut PointerMode,
        model: &mut Model,
        pid: &mut PID) {

    let pole_plot_points: egui::plot::PlotPoints = 
        model.get_poles().iter().chain(pid.get_poles().iter())
        .map(|x| x.to_cartesian()).collect();

    let zero_plot_points: egui::plot::PlotPoints = 
        model.get_zeros().iter().chain(pid.get_zeros().iter())
        .map(|x| x.to_cartesian()).collect();

    let pole_points = egui::plot::Points::new(pole_plot_points)
        .shape(egui::plot::MarkerShape::Cross)
        .radius(12.0);


    let zero_points = egui::plot::Points::new(zero_plot_points)
        .shape(egui::plot::MarkerShape::Circle)
        .radius(12.0)
        .filled(false);


    pid_config(ui, pid);
    let plot = egui::plot::Plot::new("My_plot")
        .view_aspect(1.0)
        .set_margin_fraction((0.1,0.1).into())
        .show(ui, |ui| {
        ui.points(pole_points);
        ui.points(zero_points);

        if pointer_mode == &PointerMode::Move { return }

        if ui.plot_clicked() {
            if let Some(coord) = ui.pointer_coordinate() {
                match &pointer_mode {
                    PointerMode::AddPole => model.push_pole((coord.x, coord.y)),
                    PointerMode::AddZero => model.push_zero((coord.x, coord.y)),
                    PointerMode::Remove => model.remove_closest_element((coord.x,coord.y)),
                    PointerMode::Move => {}
                }
            }
        }
    });

    let mut fill_color = [egui::Rgba::from_rgb(0.1,0.1,0.1); 4];
    let use_color = egui::Rgba::from_rgb(0.3,0.3,0.3);

    match pointer_mode {
        PointerMode::Move => { fill_color[0] = use_color },
        PointerMode::AddPole => { fill_color[1] = use_color },
        PointerMode::AddZero => { fill_color[2] = use_color },
        PointerMode::Remove => { fill_color[3] = use_color },
    }

    ui.horizontal(|ui| {
        if egui::widgets::Button::new("Move").fill(fill_color[0]).ui(ui).clicked() { 
            *pointer_mode = PointerMode::Move
        }
        if egui::widgets::Button::new("Add pole").fill(fill_color[1]).ui(ui).clicked() { 
            *pointer_mode = PointerMode::AddPole
        }
        if egui::widgets::Button::new("Add zero").fill(fill_color[2]).ui(ui).clicked() { 
            *pointer_mode = PointerMode::AddZero
        }
        if egui::widgets::Button::new("Remove").fill(fill_color[3]).ui(ui).clicked() { 
            *pointer_mode = PointerMode::Remove
        }
    });
}
