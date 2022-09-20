mod complex;
mod diagrams;
mod help;
mod pz_map;
mod response;
mod system;
mod tuner;

use diagrams::show_bode_plots;
use help::show_help;
use pz_map::show_pz_map;
use response::show_response_plot;
use tuner::show_controller_tuner;

use system::{Controller, Model};

struct OpenWindows {
    diagrams: bool,
    pz_and_tuner: bool,
    response: bool,
    help: bool,
}

impl Default for OpenWindows {
    fn default() -> Self {
        OpenWindows {
            diagrams: true,
            pz_and_tuner: true,
            help: true,
            response: true,
        }
    }
}

pub struct MyApp {
    model: Model,
    pointer_mode: pz_map::PointerMode,
    controller: Controller,
    diagram_config: diagrams::DiagramsConfiguration,
    help_config: help::HelpConfig,
    open_windows: OpenWindows,
}

impl MyApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            model: Model::default(),
            pointer_mode: pz_map::PointerMode::AddPole,
            controller: Controller::default(),
            diagram_config: diagrams::DiagramsConfiguration::default(),
            help_config: help::HelpConfig::default(),
            open_windows: OpenWindows::default(),
        }
    }
}

macro_rules! toggable {
    ($ui:expr, $variable:expr, $label:expr) => {
        if $ui.selectable_label($variable, $label).clicked() {
            $variable = !$variable
        }
    };
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let (max_x, max_y) = {
            let available_rectangle = ctx.available_rect();
            (available_rectangle.max.x, available_rectangle.max.y)
        };

        egui::containers::Window::new("Menu")
            .fixed_pos((max_x - 70.0, 0.0))
            .fixed_size((90.0, 120.0))
            .collapsible(false)
            .show(ctx, |ui| {
                toggable!(ui, self.open_windows.help, "Help");
                toggable!(ui, self.open_windows.pz_and_tuner, "Tuner");
                toggable!(ui, self.open_windows.diagrams, "Diagrams");
                toggable!(ui, self.open_windows.response, "Response");
            });

        if self.open_windows.diagrams {
            egui::containers::Window::new("Diagrams")
                .default_pos((0.0, 0.0))
                .default_size((400.0, 200.0))
                .show(ctx, |ui| {
                    show_bode_plots(ui, &self.model, &self.controller, &mut self.diagram_config);
                });
        }

        if self.open_windows.pz_and_tuner {
            egui::containers::Window::new("Tuner")
                .default_pos((max_x, 0.0))
                .show(ctx, |ui| {
                    ui.vertical(|ui| {
                        show_controller_tuner(ui, &mut self.controller);
                        show_pz_map(
                            ui,
                            &mut self.pointer_mode,
                            &mut self.model,
                            &self.controller,
                        );
                    });

                    // ui.add(egui::widgets::Label::new(format!("Poles: {:?}", self.model.get_poles())));
                    // ui.add(egui::widgets::Label::new(format!("Zeros: {:?}", self.model.get_zeros())));
                });
        }

        /*
        egui::containers::Window::new("PID transfer function")
            .default_pos((max_x / 2.0, max_y))
            .show(ctx, |ui| {
                ui.add(egui::Label::new(format!(
                    "pz: {:?}",
                    self.controller.get_pz_elements()
                )));
                ui.add(egui::Label::new(format!(
                    "a: {:?}",
                    self.controller.get_state_space().a_vector
                )));
                ui.add(egui::Label::new(format!(
                    "b: {:?}",
                    self.controller.get_state_space().b_vector
                )));
            });
        */

        if self.open_windows.response {
            egui::containers::Window::new("Response")
                .default_pos((max_x, max_y))
                .show(ctx, |ui| {
                    show_response_plot(ui, &self.model, &self.controller);
                });
        }

        if self.open_windows.help {
            egui::containers::Window::new("Help")
                .default_pos((max_x / 2.0, max_y))
                .show(ctx, |ui| {
                    show_help(ui, &mut self.help_config);
                });
        }
    }
}
