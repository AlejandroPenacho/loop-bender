use std::fmt;

pub struct HelpConfig {
    selected_menu: HelpMenu
}

impl Default for HelpConfig {
    fn default() -> Self {
        HelpConfig {
            selected_menu: HelpMenu::Main
        }
    }
}

#[derive(PartialEq,Eq)]
enum HelpMenu {
    Main,
    Diagrams,
    Tuner,
    PzMap,
    Response
}

impl fmt::Display for HelpMenu {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HelpMenu::Main => write!(f, "Main"),
            HelpMenu::Diagrams => write!(f, "Diagrams"),
            HelpMenu::Tuner => write!(f, "Tuner"),
            HelpMenu::PzMap => write!(f, "Pole-zero map"),
            HelpMenu::Response => write!(f, "Response"),
        }
    }
}


macro_rules! select_label {
    ($ui:expr, $variable:expr, $option:expr) => {
        if $ui.selectable_label($variable == $option, format!("{}", $option)).clicked() {   
            $variable = $option
        }
    }
}




pub fn show_help(ui: &mut egui::Ui, config: &mut HelpConfig) {


    ui.horizontal(|ui| {
        select_label!(ui, config.selected_menu, HelpMenu::Main);
        select_label!(ui, config.selected_menu, HelpMenu::Diagrams);
        select_label!(ui, config.selected_menu, HelpMenu::Tuner);
        select_label!(ui, config.selected_menu, HelpMenu::PzMap);
        select_label!(ui, config.selected_menu, HelpMenu::Response);
    });

    ui.separator();

    ui.label("Cuando nadie daba un duro, vengo y me cargo la escena. Jordi Wild esta en la casa, pasa que seras la cena");
}
