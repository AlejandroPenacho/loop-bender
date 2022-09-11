
pub trait DynamicalSystem {
    fn get_steady_factor(&self) -> f64;
    fn get_zeros(&self) -> Vec<(f64,f64)>;
    fn get_poles(&self) -> Vec<(f64,f64)>;
}

pub struct PID {
    k: f64,
    t_i: f64,
    t_d: f64,
    integral_on: bool,
    derivative_on: bool
}

impl PID {
    pub fn get_k(&mut self) -> &mut f64 { &mut self.k }
    pub fn get_ti(&mut self) -> &mut f64 { &mut self.t_i }
    pub fn get_td(&mut self) -> &mut f64 { &mut self.t_d }
    pub fn get_integral_on(&mut self) -> &mut bool { &mut self.integral_on }
    pub fn get_derivative_on(&mut self) -> &mut bool { &mut self.derivative_on }
}

impl DynamicalSystem for PID {
    fn get_steady_factor(&self) -> f64 { self.k }

    fn get_zeros(&self) -> Vec<(f64,f64)> {
        if !self.integral_on && !self.derivative_on {
            vec![]
        } else if self.integral_on && !self.derivative_on {
            vec![(-1.0/self.t_i, 0.0)]
        } else if !self.integral_on && self.derivative_on {
            vec![ (-1.0/self.t_d, 0.0) ]
        } else {
            let midpoint = -1.0/(2.0*self.t_d);
            let delta = 1.0/(4.0*self.t_d.powf(2.0)) - 1.0/(self.t_d * self.t_i);

            if delta >= 0.0 {
                vec![
                    (midpoint + delta.powf(0.5), 0.0),
                    (midpoint - delta.powf(0.5), 0.0),
                ]
            } else {
                vec![
                    (midpoint, (-delta).powf(0.5)),
                    (midpoint, -(-delta).powf(0.5)),
                ]
            }
        }
    }

    fn get_poles(&self) -> Vec<(f64,f64)> {
        if self.integral_on { vec![(0.0,0.0)] } else { vec![] }
    }
}

impl Default for PID {
    fn default() -> Self {
        PID { k: 1.0, t_i: 20.0, t_d: 0.0, integral_on: false, derivative_on: false }
    }
}

