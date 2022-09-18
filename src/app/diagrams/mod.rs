use crate::app::system::{Model, Controller, DynamicalSystem};

pub fn show_bode_plots(ui: &mut egui::Ui, model: &Model, controller: &Controller, bode_link: &egui::widgets::plot::LinkedAxisGroup) {

    let all_freqs_expo = (0..500).map(|i| (i as f64/100.0) - 3.0).collect::<Vec<f64>>();
    let all_freqs = all_freqs_expo.iter().map(|i| 10f64.powf(*i)).collect::<Vec<f64>>();

    let mut phase_cross_point = None;
    let mut gain_cross_point = None;


    let mut prev_mag = None;

    let mags = all_freqs.iter().map(|f| {
        let mut mag = 1f64;

        mag *= model.get_mag_at_freq(*f);
        mag *= controller.get_mag_at_freq(*f);

        let log_mag = mag.log10();

        if prev_mag.map_or(false, |prev| prev > 0.0 && log_mag <= 0.0 ) && gain_cross_point.is_none() {
            gain_cross_point = Some(f);
        }
        prev_mag = Some(log_mag);

        log_mag
    }).collect::<Vec<f64>>();
    
    let mut prev_phase = None;
    let phases = all_freqs.iter().map(|f| {
        let mut phase = 0.0;

        phase += model.get_phase_at_freq(*f);
        phase += controller.get_phase_at_freq(*f);

        if prev_phase.map_or(false, |prev| prev >= -std::f64::consts::PI && phase < -std::f64::consts::PI) && phase_cross_point.is_none() {
            phase_cross_point = Some(f);
        }
        prev_phase = Some(phase);

        phase
    }).collect::<Vec<f64>>();

    let phase_margin = gain_cross_point.map(|x| std::f64::consts::PI + (model.get_phase_at_freq(*x) + controller.get_phase_at_freq(*x)));
    let gain_margin = phase_cross_point.map(|x| -(model.get_mag_at_freq(*x)*controller.get_mag_at_freq(*x)).log10());

    let mag_line: egui::plot::PlotPoints = (0..500).map(|i| [all_freqs_expo[i], mags[i]]).collect();
    let phase_line: egui::plot::PlotPoints  = (0..500).map(|i| [all_freqs_expo[i], phases[i]]).collect();

    let nyquist_line: egui::plot::PlotPoints =
        (0..500).map(|i| [10f64.powf(mags[i])*phases[i].cos(), 10f64.powf(mags[i])*phases[i].sin()])
        .collect();

    let phase_margin_label = match phase_margin {
        Some(x) => format!("{:.1}", x*180.0/std::f64::consts::PI),
        None => "∞".to_owned()
    };

    let gain_margin_label = match gain_margin {
        Some(x) => format!("{:.1}", 10.0*x),
        None => "∞".to_owned()
    };

    ui.vertical(|ui| {
        ui.add(egui::widgets::Label::new(format!("Gain margin: {} dB", gain_margin_label)));
        ui.add(egui::widgets::Label::new(format!("Phase margin: {} °", phase_margin_label)));
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
                            ui.vline(egui::plot::VLine::new(freq.log10()));
                        }
                        if let Some(freq) = phase_cross_point {
                            ui.vline(egui::plot::VLine::new(freq.log10()));
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
                            ui.vline(egui::plot::VLine::new(freq.log10()));
                        }
                        if let Some(freq) = phase_cross_point {
                            ui.vline(egui::plot::VLine::new(freq.log10()));
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
