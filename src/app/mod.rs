mod complex;
mod system;
mod tuner;
mod pz_map;
mod diagrams;
mod response;

use tuner::show_controller_tuner;
use pz_map::show_pz_map;
use diagrams::show_bode_plots;
use response::show_response_plot;

use system::{Controller, Model, DynamicalSystem};

pub struct MyApp {
    model: Model,
    pointer_mode: pz_map::PointerMode,
    controller: Controller,
    diagram_config: diagrams::DiagramsConfiguration
}

impl MyApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            model: Model::default(),
            pointer_mode: pz_map::PointerMode::AddPole,
            controller: Controller::default(),
            diagram_config: diagrams::DiagramsConfiguration::default()
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

            show_bode_plots(ui, &self.model, &self.controller, &mut self.diagram_config);
        });

        egui::containers::Window::new("Tuner")
            .default_pos((max_x,0.0))
            .show(ctx, |ui| {

            ui.vertical(|ui| {
                show_controller_tuner(ui, &mut self.controller);
                show_pz_map(ui, &mut self.pointer_mode, &mut self.model, &self.controller);
            });

            // ui.add(egui::widgets::Label::new(format!("Poles: {:?}", self.model.get_poles())));
            // ui.add(egui::widgets::Label::new(format!("Zeros: {:?}", self.model.get_zeros())));
        });

        egui::containers::Window::new("PID transfer function")
            .default_pos((max_x/2.0, max_y))
            .show(ctx, |ui| {
                ui.add(egui::Label::new(format!("pz: {:?}", self.controller.get_pz_elements())));
                ui.add(egui::Label::new(format!("a: {:?}", self.controller.get_state_space().a_vector)));
                ui.add(egui::Label::new(format!("b: {:?}", self.controller.get_state_space().b_vector)));
            });

        egui::containers::Window::new("Response")
            .default_pos((max_x, max_y))
            .show(ctx, |ui| {
                show_response_plot(ui, &self.model, &self.controller);
            });
    }
}

