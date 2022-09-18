use crate::app::system::{Model, Controller, DynamicalSystem};

pub struct DiagramsConfiguration {
    show_gain_margin: bool,
    show_phase_margin: bool,
    bode_link: egui::widgets::plot::LinkedAxisGroup
}

impl Default for DiagramsConfiguration {
    fn default() -> Self {
        DiagramsConfiguration {
            show_phase_margin: false,
            show_gain_margin: false,
            bode_link: egui::widgets::plot::LinkedAxisGroup::new(true, false)
        }
    }
}

struct StabilityMargins {
    gain_margin: Option<(f64, f64)>,
    phase_margin: Option<(f64, f64)>,
}

pub fn show_bode_plots(
    ui: &mut egui::Ui,
    model: &Model,
    controller: &Controller,
    config: &mut DiagramsConfiguration) {

    let all_freqs_expo = (0..500).map(|i| (i as f64/100.0) - 3.0).collect::<Vec<f64>>();
    let all_freqs = all_freqs_expo.iter().map(|i| 10f64.powf(*i)).collect::<Vec<f64>>();
    
    let system = model.link_system(controller);

    let (mags, phases, margins) = compute_phase_and_margin(&system, &all_freqs);

    let mag_line: egui::plot::PlotPoints = (0..500).map(|i| [all_freqs_expo[i], mags[i]]).collect();
    let phase_line: egui::plot::PlotPoints  = (0..500).map(|i| [all_freqs_expo[i], phases[i]]).collect();

    let nyquist_line: egui::plot::PlotPoints =
        (0..500).map(|i| [10f64.powf(mags[i])*phases[i].cos(), 10f64.powf(mags[i])*phases[i].sin()])
        .collect();

    let phase_margin_label = match margins.phase_margin {
        Some(x) => format!("{:.1}", x.1*180.0/std::f64::consts::PI),
        None => "∞".to_owned()
    };

    let gain_margin_label = match margins.gain_margin {
        Some(x) => format!("{:.1}", 10.0*x.1),
        None => "∞".to_owned()
    };

    ui.vertical(|ui| {
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.add(egui::widgets::Label::new("Gain margin:"));
                ui.add(egui::widgets::Label::new("Phase margin:"));
            });
            ui.vertical(|ui| {
                ui.add(egui::widgets::Label::new(format!("{:<6} dB", gain_margin_label)));
                ui.add(egui::widgets::Label::new(format!("{:<6} °", phase_margin_label)));
            });
            ui.vertical(|ui| {
                ui.add(egui::widgets::Checkbox::new(&mut config.show_gain_margin, "Show"));
                ui.add(egui::widgets::Checkbox::new(&mut config.show_phase_margin, "Show"));
            });
        });
        ui.horizontal(|ui| {
            ui.vertical(|ui| {
                ui.add(egui::widgets::Label::new("Magnitude"));
                egui::plot::Plot::new("bode_mag")
                    .view_aspect(2.0)
                    .center_y_axis(false)
                    .width(ui.available_width()/2.1)
                    .link_axis(config.bode_link.clone())
                    .x_axis_formatter(|x,_| if x < 3.5 && x > -3.5 {format!("{}", 10f64.powf(x))} else {format!("{:e}", 10f64.powf(x))})
                    .y_axis_formatter(|x,_| format!("{} dB", x*10.0))
                    .show(ui,
                    |ui| {
                        ui.line(egui::plot::Line::new(mag_line));
                        if let Some(freq) = margins.gain_margin {
                            if config.show_gain_margin {
                                ui.vline(egui::plot::VLine::new(freq.0.log10()));
                            }
                        }
                        if let Some(freq) = margins.phase_margin {
                            if config.show_phase_margin {
                                ui.vline(egui::plot::VLine::new(freq.0.log10()));
                            }
                        }
                    });

                ui.add(egui::widgets::Label::new("Phase"));
                egui::plot::Plot::new("bode_phase")
                    .view_aspect(2.0)
                    .center_y_axis(false)
                    .width(ui.available_width()/2.1)
                    .link_axis(config.bode_link.clone())
                    .x_axis_formatter(|x,_| if x < 3.5 && x > -3.5 {format!("{}", 10f64.powf(x))} else {format!("{:e}", 10f64.powf(x))})
                    .y_grid_spacer(bode_phase_y_spacer)
                    .y_axis_formatter(bode_phase_y_formatter)
                    .show(ui,
                    |ui| {
                        ui.line(egui::plot::Line::new(phase_line));
                        if let Some(freq) = margins.gain_margin {
                            if config.show_gain_margin {
                                ui.vline(egui::plot::VLine::new(freq.0.log10()));
                            }
                        }
                        if let Some(freq) = margins.phase_margin {
                            if config.show_phase_margin {
                                ui.vline(egui::plot::VLine::new(freq.0.log10()));
                            }
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

fn compute_phase_and_margin<T: DynamicalSystem>(
    system: &T,
    frequencies: &[f64]) -> (Vec<f64>, Vec<f64>, StabilityMargins) {

    let mut phase_cross_point = None;
    let mut gain_cross_point = None;


    let mut prev_mag = None;

    let mags = frequencies.iter().map(|f| {
        let mut mag = 1f64;

        mag *= system.get_mag_at_freq(*f);

        let log_mag = mag.log10();

        if prev_mag.map_or(false, |prev| prev > 0.0 && log_mag <= 0.0 ) && gain_cross_point.is_none() {
            gain_cross_point = Some(f);
        }
        prev_mag = Some(log_mag);

        log_mag
    }).collect::<Vec<f64>>();
    
    let mut prev_phase = None;
    let phases = frequencies.iter().map(|f| {
        let mut phase = 0.0;

        phase += system.get_phase_at_freq(*f);

        if prev_phase.map_or(false, |prev| prev >= -std::f64::consts::PI && phase < -std::f64::consts::PI) && phase_cross_point.is_none() {
            phase_cross_point = Some(f);
        }
        prev_phase = Some(phase);

        phase
    }).collect::<Vec<f64>>();

    let margins = StabilityMargins {
        gain_margin: phase_cross_point.map(|x| (*x, -system.get_mag_at_freq(*x).log10())),
        phase_margin: gain_cross_point.map(|x| (*x, std::f64::consts::PI + (system.get_phase_at_freq(*x))))
    };

    (mags, phases, margins)

}
