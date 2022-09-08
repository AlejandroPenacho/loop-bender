
fn main() {
    eframe::run_native("My first App", eframe::NativeOptions::default(), Box::new(|cc| Box::new(MyApp::new(cc))));
}

struct MyApp {
    counter: i32,
    poles: Vec<(f64,f64)>,
    zeros: Vec<(f64,f64)>,
    pointer_mode: PointerMode
}

enum PointerMode {
    AddPole,
    AddZero,
    Remove
}

impl PointerMode {
    fn display(&self) -> &'static str {
        match self {
            PointerMode::AddPole => "Pole",
            PointerMode::AddZero => "Zero",
            PointerMode::Remove => "Remove"
        }
    }
}


impl MyApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            poles: vec![(1.0,1.0), (2.0,1.0)],
            zeros: vec![(-2.0, 0.5)],
            counter: 0,
            pointer_mode: PointerMode::AddPole
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Heyyy");
            ui.vertical(|ui| {
                ui.horizontal(|ui| {
                    if ui.button(self.pointer_mode.display()).clicked() {
                        match self.pointer_mode {
                            PointerMode::AddPole => self.pointer_mode = PointerMode::AddZero,
                            PointerMode::AddZero => self.pointer_mode = PointerMode::Remove,
                            PointerMode::Remove => self.pointer_mode = PointerMode::AddPole
                        }
                    }
                });
                get_pz_map(ui, &self.pointer_mode, &mut self.poles, &mut self.zeros);
            })
        });
    }
}

fn get_pz_map(ui: &mut egui::Ui, pointer_mode: &PointerMode, poles: &mut Vec<(f64, f64)>, zeros: &mut Vec<(f64, f64)>) {

    let pole_plot_points: egui::plot::PlotPoints = poles.iter().map(|(x,y)| [*x, *y]).collect();
    let zero_plot_points: egui::plot::PlotPoints = zeros.iter().map(|(x,y)| [*x, *y]).collect();

    let pole_points = egui::plot::Points::new(pole_plot_points).shape(egui::plot::MarkerShape::Cross).radius(12.0);
    let zero_points = egui::plot::Points::new(zero_plot_points).shape(egui::plot::MarkerShape::Circle).radius(12.0);

    let plot = egui::plot::Plot::new("My_plot").view_aspect(1.0).show(ui,
        |ui| {
            ui.points(pole_points);
            ui.points(zero_points);

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

                        }
                    }

                }
            }
    });

}
