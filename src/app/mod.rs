mod complex;
mod system;

use system::{PID, Model, DynamicalSystem};

use egui::widgets::Widget;

#[derive(Default)]
struct StabilityMargins {
    gain_margin: Option<f64>,
    phase_margin: Option<f64>
}

pub struct MyApp {
    model: Model,
    pointer_mode: PointerMode,
    pid: PID,
    bode_link: egui::widgets::plot::LinkedAxisGroup,
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
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            model: Model::default(),
            pointer_mode: PointerMode::AddPole,
            pid: PID::default(),
            bode_link: egui::widgets::plot::LinkedAxisGroup::new(true, false),
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let (max_x, max_y) = {
            let available_rectangle = ctx.available_rect();
            (available_rectangle.max.x, available_rectangle.max.y)
        };

        egui::containers::Window::new("Diagrams")
            .default_pos((0.0,0.0))
            .default_size((400.0,200.0))
            .show(ctx, |ui| {

            get_bode_plot(ui, &self.model, &self.pid, &self.bode_link);
        });

        egui::containers::Window::new("Tuner")
            .default_pos((max_x,0.0))
            .show(ctx, |ui| {

            get_pz_map(ui, &mut self.pointer_mode, &mut self.model, &mut self.pid);
            // ui.add(egui::widgets::Label::new(format!("Poles: {:?}", self.model.get_poles())));
            // ui.add(egui::widgets::Label::new(format!("Zeros: {:?}", self.model.get_zeros())));
        });

        egui::containers::Window::new("PID transfer function")
            .default_pos((max_x/2.0, max_y))
            .show(ctx, |ui| {
                ui.add(egui::Label::new(format!("pz: {:?}", self.pid.get_pz_elements())));
                ui.add(egui::Label::new(format!("a: {:?}", self.pid.get_state_space().a_vector)));
                ui.add(egui::Label::new(format!("b: {:?}", self.pid.get_state_space().b_vector)));
            });

        egui::containers::Window::new("Response")
            .default_pos((max_x, max_y))
            .show(ctx, |ui| {
                response_plot(ui, &self.model, &self.pid);
            });
    }
}

fn response_plot(ui: &mut egui::Ui, model: &Model, pid: &PID) {
    let mut state_space_model = model.get_state_space();
    let mut state_space_pid = pid.get_state_space();

    let mut output: Vec<f64> = Vec::with_capacity(5000);
    let time = (0..5000).map(|x| x as f64/500.0).collect::<Vec<f64>>();

    output.push(0.0);

    for i in 0..5000 {
        /*
        state_space_pid.step(0.01, i as f64 / 1000.0);
        output.push(state_space_pid.get_output());
        */

        let model_output = state_space_model.get_output();
        let pid_output = state_space_pid.get_output();

        state_space_model.step(0.002, pid_output);
        state_space_pid.step(0.002, 1.0-model_output);
        output.push(state_space_model.get_output());
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

    let mut phase_cross_point = None;
    let mut gain_cross_point = None;
    let mut prev_mag = None;

    let mags = all_freqs.iter().map(|f| {
        let mut mag = 1f64;

        mag *= model.get_mag_at_freq(*f);
        mag *= pid.get_mag_at_freq(*f);

        let log_mag = mag.log10();

        if prev_mag.map_or(false, |prev| prev > 0.0 && log_mag <= 0.0 ) && gain_cross_point.is_none() {
            gain_cross_point = Some(f.log10());
        }
        prev_mag = Some(log_mag);

        log_mag
    }).collect::<Vec<f64>>();
    
    let mut prev_phase = None;
    let phases = all_freqs.iter().map(|f| {
        let mut phase = 0.0;

        phase += model.get_phase_at_freq(*f);
        phase += pid.get_phase_at_freq(*f);

        if prev_phase.map_or(false, |prev| prev >= -std::f64::consts::PI && phase < -std::f64::consts::PI) && phase_cross_point.is_none() {
            phase_cross_point = Some(f.log10());
        }
        prev_phase = Some(phase);

        phase
    }).collect::<Vec<f64>>();

    let mag_line: egui::plot::PlotPoints = (0..500).map(|i| [all_freqs_expo[i], mags[i]]).collect();
    let phase_line: egui::plot::PlotPoints  = (0..500).map(|i| [all_freqs_expo[i], phases[i]]).collect();

    let nyquist_line: egui::plot::PlotPoints =
        (0..500).map(|i| [10f64.powf(mags[i])*phases[i].cos(), 10f64.powf(mags[i])*phases[i].sin()])
        .collect();


    ui.horizontal(|ui| {
        ui.vertical(|ui| {
            ui.add(egui::widgets::Label::new("Magnitude"));
            egui::plot::Plot::new("bode_mag")
                .view_aspect(2.0)
                .center_y_axis(false)
                .width(ui.available_width()/2.1)
                .link_axis(bode_link.clone())
                .x_axis_formatter(|x,_| if x < 3.5 && x > -3.5 {format!("{}", 10f64.powf(x))} else {format!("{:e}", 10f64.powf(x))})
                .y_axis_formatter(|x,_| format!("{} dB", x*10.0))
                .show(ui,
                |ui| {
                    ui.line(egui::plot::Line::new(mag_line));
                    if let Some(freq) = gain_cross_point {
                        ui.vline(egui::plot::VLine::new(freq));
                    }
                    if let Some(freq) = phase_cross_point {
                        ui.vline(egui::plot::VLine::new(freq));
                    }
                });

            ui.add(egui::widgets::Label::new("Phase"));
            egui::plot::Plot::new("bode_phase")
                .view_aspect(2.0)
                .center_y_axis(false)
                .width(ui.available_width()/2.1)
                .link_axis(bode_link.clone())
                .x_axis_formatter(|x,_| if x < 3.5 && x > -3.5 {format!("{}", 10f64.powf(x))} else {format!("{:e}", 10f64.powf(x))})
                .y_grid_spacer(bode_phase_y_spacer)
                .y_axis_formatter(bode_phase_y_formatter)
                .show(ui,
                |ui| {
                    ui.line(egui::plot::Line::new(phase_line));
                    if let Some(freq) = gain_cross_point {
                        ui.vline(egui::plot::VLine::new(freq));
                    }
                    if let Some(freq) = phase_cross_point {
                        ui.vline(egui::plot::VLine::new(freq));
                    }
                });
        });

        ui.vertical(|ui| {
            ui.add(egui::widgets::Label::new("Nyquist"));
            egui::plot::Plot::new("nyquist")
                .view_aspect(1.0)
                .width(ui.available_width())
                .data_aspect(1.0)
                .include_x(0.0)
                .include_y(0.0)
                .show(ui,
                |ui| {
                    ui.line(egui::plot::Line::new(nyquist_line));
                }
            );
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

fn bode_phase_y_spacer(grid: egui::plot::GridInput) -> Vec<egui::plot::GridMark> {
    let mut output = Vec::new();
    for i in -8..8 {
        output.push(egui::plot::GridMark {
            value: i as f64/4.0 * std::f64::consts::PI,
            step_size:  1.0/4.0 * std::f64::consts::PI,
        });
    }

    output
}

fn bode_phase_y_formatter(y: f64, _range: &std::ops::RangeInclusive<f64>) -> String {
    let quarters = (4.0 * y / std::f64::consts::PI).round() as i32;

    if quarters % 4 == 0 {
        format!("{} π", quarters/4)
    } else if quarters % 2 == 0 {
        format!("{}/2 π", quarters/2)
    } else {
        format!("{}/4 π", quarters)
    }
}
