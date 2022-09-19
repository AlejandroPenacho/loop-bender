use crate::app::system::Controller;

pub fn show_controller_tuner(ui: &mut egui::Ui, controller: &mut Controller) {
    let pid = &mut controller.pid;
    let lead_lag = &mut controller.lead_lag;
    let lag_lead = &mut controller.lag_lead;

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

        let k_speed = *pid.get_k() / 20.0;
        let i_speed = *pid.get_ti() / 20.0;
        let d_speed = *pid.get_td() / 20.0;

        ui.vertical(|ui| {
            ui.add(
                egui::widgets::DragValue::new(pid.get_k())
                    .clamp_range(0.01..=50.0)
                    .speed(k_speed),
            );

            ui.add_enabled(
                *pid.get_integral_on(),
                egui::widgets::DragValue::new(pid.get_ti())
                    .clamp_range(0.1..=50.0)
                    .speed(i_speed),
            );

            ui.add_enabled(
                *pid.get_derivative_on(),
                egui::widgets::DragValue::new(pid.get_td())
                    .clamp_range(0.1..=50.0)
                    .speed(d_speed),
            );
        });

        ui.add(egui::widgets::Separator::default().spacing(30.0));

        let lead_lag_zero_speed = (*lead_lag.get_zero_freq() / 20.0).max(0.0001);
        let lead_lag_pole_speed = (*lead_lag.get_pole_freq() / 20.0).max(0.0001);
        let pole_value = *lead_lag.get_pole_freq();
        let zero_value = *lead_lag.get_zero_freq();

        ui.vertical(|ui| {
            ui.add(egui::widgets::Checkbox::new(
                lead_lag.is_activated(),
                "Lead-lag",
            ));
            ui.horizontal(|ui| {
                ui.add_enabled(*lead_lag.is_activated(), egui::widgets::Label::new("Zero:"));
                ui.add_enabled(
                    *lead_lag.is_activated(),
                    egui::widgets::DragValue::new(lead_lag.get_zero_freq())
                        .clamp_range(0.0..=pole_value)
                        .speed(lead_lag_zero_speed),
                );
            });
            ui.horizontal(|ui| {
                ui.add_enabled(*lead_lag.is_activated(), egui::widgets::Label::new("Pole:"));
                ui.add_enabled(
                    *lead_lag.is_activated(),
                    egui::widgets::DragValue::new(lead_lag.get_pole_freq())
                        .clamp_range(zero_value..=10000.0)
                        .speed(lead_lag_pole_speed),
                );
            });
        });

        ui.add(egui::widgets::Separator::default().spacing(30.0));

        let lag_lead_zero_speed = (*lag_lead.get_zero_freq() / 20.0).max(0.0001);
        let lag_lead_pole_speed = (*lag_lead.get_pole_freq() / 20.0).max(0.0001);
        let pole_value = *lag_lead.get_pole_freq();
        let zero_value = *lag_lead.get_zero_freq();

        ui.vertical(|ui| {
            ui.add(egui::widgets::Checkbox::new(
                lag_lead.is_activated(),
                "Lag-lead",
            ));
            ui.horizontal(|ui| {
                ui.add_enabled(*lag_lead.is_activated(), egui::widgets::Label::new("Zero:"));
                ui.add_enabled(
                    *lag_lead.is_activated(),
                    egui::widgets::DragValue::new(lag_lead.get_zero_freq())
                        .clamp_range(pole_value..=10000.0)
                        .speed(lag_lead_zero_speed),
                );
            });
            ui.horizontal(|ui| {
                ui.add_enabled(*lag_lead.is_activated(), egui::widgets::Label::new("Pole:"));
                ui.add_enabled(
                    *lag_lead.is_activated(),
                    egui::widgets::DragValue::new(lag_lead.get_pole_freq())
                        .clamp_range(0.0..=zero_value)
                        .speed(lag_lead_pole_speed),
                );
            });
        });
    });
}
