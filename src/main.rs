mod complex;
mod system;

use system::{PID, DynamicalSystem};

use egui::widgets::Widget;

fn main() {
    eframe::run_native("My first App", eframe::NativeOptions::default(), Box::new(|cc| Box::new(MyApp::new(cc))));
}


struct MyApp {
    poles: Vec<(f64,f64)>,
    zeros: Vec<(f64,f64)>,
    pointer_mode: PointerMode,
    pid: PID
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
            poles: vec![],
            zeros: vec![],
            pointer_mode: PointerMode::AddPole,
            pid: PID::default()
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::containers::Window::new("left")
            .default_pos((0.0,0.0))
            .show(ctx, |ui| {

            ui.heading("Heyyy");
            get_bode_plot(ui, &self.poles, &self.zeros, &mut self.pid);
        });

        egui::containers::Window::new("right")
            .default_pos((100.0,0.0))
            .show(ctx, |ui| {

            get_pz_map(ui, &mut self.pointer_mode, &mut self.poles, &mut self.zeros, &mut self.pid);
        });
    }
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

        let k_speed = *pid.get_k()/5.0;
        let i_speed = *pid.get_ti()/5.0;
        let d_speed = *pid.get_td()/5.0;

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

fn get_bode_plot(ui: &mut egui::Ui, poles: &Vec<(f64, f64)>, zeros: &Vec<(f64,f64)>, pid: &mut PID) {
    let all_freqs_expo = (0..500).map(|i| (i as f64/100.0) - 3.0).collect::<Vec<f64>>();
    let all_freqs = all_freqs_expo.iter().map(|i| 10f64.powf(*i)).collect::<Vec<f64>>();

    let mags = all_freqs.iter().map(|f| {
        let mut mag = 0f64;
        for pole in poles.iter().chain(pid.get_poles().iter()) {
            mag -= 0.5 * (pole.0.powf(2.0) + (pole.1 - f).powf(2.0)).log10();
            if pole.0 != 0.0 && pole.0 != 0.0 {
                mag += 0.5 * (pole.0.powf(2.0) + pole.1.powf(2.0)).log10();
            }
        }

        for zero in zeros.iter().chain(pid.get_zeros().iter()) {
            mag += 0.5 * (zero.0.powf(2.0) + (zero.1 - f).powf(2.0)).log10();
            if zero.0 != 0.0 && zero.0 != 0.0 {
                mag -=  0.5 * (zero.0.powf(2.0) + zero.1.powf(2.0)).log10();
            }
        }
        mag + pid.get_steady_factor().log10()
    }).collect::<Vec<f64>>();
    
    let phases = all_freqs.iter().map(|f| {
        let mut phase = 0f64;
        for pole in poles.iter() .chain(pid.get_poles().iter()) {
            let delta_phase = (f - pole.1).atan2(0.0 - pole.0);
            /*
            if delta_phase < 0.0 {
                delta_phase = 2.0*std::f64::consts::PI + delta_phase;
            }
            */
            phase -= delta_phase;
        }
        for zero in zeros.iter().chain(pid.get_zeros().iter()) {
            let delta_phase = (f - zero.1).atan2(0.0 - zero.0);
            /*
            if delta_phase < 0.0 {
                delta_phase = 2.0*std::f64::consts::PI + delta_phase;
            }
            */
            phase += delta_phase;
        }
        phase
    }).collect::<Vec<f64>>();

    let mag_line: egui::plot::PlotPoints = (0..500).map(|i| [all_freqs_expo[i], mags[i]]).collect();
    let phase_line: egui::plot::PlotPoints  = (0..500).map(|i| [all_freqs_expo[i], phases[i]]).collect();

    let nyquist_line: egui::plot::PlotPoints = (0..500).map(|i| [10f64.powf(mags[i])*phases[i].cos(), 10f64.powf(mags[i])*phases[i].sin()]).collect();

    ui.vertical(|ui| {
        // ui.add(egui::widgets::Label::new(format!("{:?}", mags)));
        egui::plot::Plot::new("bode_mag")
            .view_aspect(2.0)
            .center_y_axis(false)
            .show(ui,
            |ui| {
                ui.line(egui::plot::Line::new(mag_line));
            });
        egui::plot::Plot::new("bode_phase")
            .view_aspect(2.0)
            .center_y_axis(false)
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
        poles: &mut Vec<(f64, f64)>,
        zeros: &mut Vec<(f64, f64)>,
        pid: &mut PID) {

    let pole_plot_points: egui::plot::PlotPoints = poles.iter().chain(pid.get_poles().iter()).map(|(x,y)| [*x, *y]).collect();
    let zero_plot_points: egui::plot::PlotPoints = zeros.iter().chain(pid.get_zeros().iter()).map(|(x,y)| [*x, *y]).collect();

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
                match pointer_mode {
                    PointerMode::AddPole => poles.push((coord.x, coord.y)),
                    PointerMode::AddZero => zeros.push((coord.x, coord.y)),
                    PointerMode::Remove => {
                        let closest_pole = poles.iter()
                            .map(|(x,y)| (x - coord.x).powf(2.0) + (y - coord.y).powf(2.0))
                            .enumerate().reduce(|(i_min, dist_min), (i_new, dist_new)| {
                                if let Some(std::cmp::Ordering::Less) = dist_new.partial_cmp(&dist_min) {
                                    (i_new, dist_new)
                                } else {
                                    (i_min, dist_min)
                                }
                            });
                        let closest_zero = zeros.iter()
                            .map(|(x,y)| (x - coord.x).powf(2.0) + (y - coord.y).powf(2.0))
                            .enumerate().reduce(|(i_min, dist_min), (i_new, dist_new)| {
                                if let Some(std::cmp::Ordering::Less) = dist_new.partial_cmp(&dist_min) {
                                    (i_new, dist_new)
                                } else {
                                    (i_min, dist_min)
                                }
                            });

                        if closest_zero.is_none() && closest_pole.is_none() { };
                        if closest_zero.is_none() && closest_pole.is_some() { poles.remove(closest_pole.unwrap().0); };
                        if closest_zero.is_some() && closest_pole.is_none() { zeros.remove(closest_zero.unwrap().0); };
                        if closest_zero.is_some() && closest_pole.is_some() {
                            if closest_pole.unwrap().1 < closest_zero.unwrap().1 {
                                poles.remove(closest_pole.unwrap().0) ;
                            } else {
                                zeros.remove(closest_zero.unwrap().0);
                            }
                        };
                    },
                    _ => {}
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
