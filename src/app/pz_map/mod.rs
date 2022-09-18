use crate::app::system::{Model, Controller, DynamicalSystem};

use egui::widgets::Widget;

#[derive(PartialEq,Eq)]
pub enum PointerMode {
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

pub fn show_pz_map(
        ui: &mut egui::Ui,
        pointer_mode: &mut PointerMode,
        model: &mut Model,
        controller: &Controller) {

    let pid = &controller.pid;

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
