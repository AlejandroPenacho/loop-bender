mod app;

use app::MyApp;

#[cfg(not(target_arch = "wasm32"))]
fn main() {
    eframe::run_native(
        "My first App",
        eframe::NativeOptions::default(),
        Box::new(|cc| Box::new(MyApp::new(cc))),
    );
}

#[cfg(target_arch = "wasm32")]
fn main() {
    console_error_panic_hook::set_once();
    tracing_wasm::set_as_global_default();

    eframe::start_web(
        "the_canvas_id",
        eframe::WebOptions::default(),
        Box::new(|cc| Box::new(MyApp::new(cc))),
    )
    .expect("Failed");
}
