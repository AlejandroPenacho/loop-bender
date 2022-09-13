use super::complex::ComplexNumber;

pub trait DynamicalSystem {
    fn get_steady_factor(&self) -> f64;
    fn get_zeros(&self) -> Vec<ComplexNumber>;
    fn get_poles(&self) -> Vec<ComplexNumber>;

    fn get_mag_at_freq(&self, freq: f64) -> f64 {
        let complex_freq = ComplexNumber::from_cartesian(0.0, freq);
        let mut magnitude = 1.0;
        for zero in self.get_zeros() {
            magnitude *= (complex_freq - zero).get_mag();
            if !zero.is_origin() {
                magnitude /= zero.get_mag();
            }
        }
        for pole in self.get_poles() {
            magnitude /= (complex_freq - pole).get_mag();
            if !pole.is_origin() {
                magnitude *= pole.get_mag();
            }
        }
        magnitude * self.get_steady_factor()
    }

    fn get_phase_at_freq(&self, freq: f64) -> f64 {
        let complex_freq = ComplexNumber::from_cartesian(0.0, freq);
        let mut phase = 0.0;
        for zero in self.get_zeros() {
            phase += (complex_freq - zero).get_phase();
        }
        for pole in self.get_poles() {
            phase -= (complex_freq - pole).get_phase();
        }
        phase
    }

    fn get_state_space(&self) -> StateSpace {
        StateSpace {
            n_states: 1,
            state_values: vec![0.0],
            a_vector: vec![1.0],
            b_vector: vec![1.0, 1.0]
        }
    }
}

impl DynamicalSystem for Vec<Box<dyn DynamicalSystem>> {
    fn get_steady_factor(&self) -> f64 {
        self.iter().map(|x| x.get_steady_factor()).product::<f64>()
    }
    fn get_zeros(&self) -> Vec<ComplexNumber> {
        self.iter().map(|x| x.get_zeros()).flatten().collect()
    }
    fn get_poles(&self) -> Vec<ComplexNumber> {
        self.iter().map(|x| x.get_poles()).flatten().collect()
    }
}

#[derive(Debug,Clone)]
enum PZElement {
    Pole(ComplexNumber),
    DoublePole(ComplexNumber),
    Zero(ComplexNumber),
    DoubleZero(ComplexNumber)
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

    fn get_zeros(&self) -> Vec<ComplexNumber> {
        if !self.integral_on && !self.derivative_on {
            vec![]
        } else if self.integral_on && !self.derivative_on {
            vec![ComplexNumber::from_cartesian(-1.0/self.t_i, 0.0)]
        } else if !self.integral_on && self.derivative_on {
            vec![ ComplexNumber::from_cartesian(-1.0/self.t_d, 0.0) ]
        } else {
            let midpoint = -1.0/(2.0*self.t_d);
            let delta = 1.0/(4.0*self.t_d.powf(2.0)) - 1.0/(self.t_d * self.t_i);

            if delta >= 0.0 {
                vec![
                    ComplexNumber::from_cartesian(midpoint + delta.powf(0.5), 0.0),
                    ComplexNumber::from_cartesian(midpoint - delta.powf(0.5), 0.0),
                ]
            } else {
                vec![
                    ComplexNumber::from_cartesian(midpoint, (-delta).powf(0.5)),
                    ComplexNumber::from_cartesian(midpoint, -(-delta).powf(0.5)),
                ]
            }
        }
    }

    fn get_poles(&self) -> Vec<ComplexNumber> {
        if self.integral_on { vec![ComplexNumber::from_cartesian(0.0,0.0)] } else { vec![] }
    }
}

impl Default for PID {
    fn default() -> Self {
        PID { k: 1.0, t_i: 20.0, t_d: 0.0, integral_on: false, derivative_on: false }
    }
}

#[derive(Default)]
pub struct Model {
    components: Vec<PZElement>
}

impl Model {
    pub fn push_pole(&mut self, coord: (f64, f64)) {
        if coord.1.abs() < coord.0.abs()/10.0 {
            self.components.push(PZElement::Pole(ComplexNumber::from_cartesian(coord.0, 0.0)));
        } else {
            self.components.push(PZElement::DoublePole(ComplexNumber::from_cartesian(coord.0, coord.1)));
        }
    }

    pub fn push_zero(&mut self, coord: (f64, f64)) {
        if coord.1.abs() < coord.0.abs()/10.0 {
            self.components.push(PZElement::Zero(ComplexNumber::from_cartesian(coord.0, 0.0)));
        } else {
            self.components.push(PZElement::DoubleZero(ComplexNumber::from_cartesian(coord.0, coord.1)));
        }
    }
    pub fn remove_closest_element(&mut self, coord: (f64, f64)) {
        let x = coord.0;
        let y = coord.1;

        let closest_element = self.components.iter().map(|element| {
            let distance = match element {
                PZElement::Zero(point) | PZElement::Pole(point) => {
                    let [p_x, p_y] = point.to_cartesian();
                    (p_x-x).powf(2.0) + (p_y-y).powf(2.0)
                },
                PZElement::DoubleZero(point) | PZElement::DoublePole(point) => {
                    let [p_x, p_y] = point.to_cartesian();
                    (p_x-x).powf(2.0) + (p_y-y).powf(2.0).min((p_y+y).powf(2.0))

                }
            };
            distance
        }).enumerate().reduce(|(i_min, dist_min), (i_new, dist_new)| {
            if let Some(std::cmp::Ordering::Less) = dist_new.partial_cmp(&dist_min) {
                (i_new, dist_new)
            } else {
                (i_min, dist_min)
            }
        }).map(|x| x.0);

        if let Some(index) = closest_element {
            self.components.remove(index);
        }
    }
}

impl DynamicalSystem for Model {
    fn get_steady_factor(&self) -> f64 { 1.0 }

    fn get_poles(&self) -> Vec<ComplexNumber> {
        self.components.iter().filter_map(|x| {
            match x {
                PZElement::Pole(p) => Some(vec![*p].into_iter()),
                PZElement::DoublePole(p) => {
                    let conjugate_pole = p.conjugate();
                    Some(vec![*p, conjugate_pole].into_iter())
                },
                _ => None
            }
        }).flatten().collect()
    }
    fn get_zeros(&self) -> Vec<ComplexNumber> {
        self.components.iter().filter_map(|x| {
            match x {
                PZElement::Zero(z) => Some(vec![*z].into_iter()),
                PZElement::DoubleZero(z) => {
                    let conjugate_zero = z.conjugate();
                    Some(vec![*z, conjugate_zero].into_iter())
                },
                _ => None
            }
        }).flatten().collect()
    }
}

/// Represents a set of differential equations that allow the simulation of a
/// dynamic system
///
/// The transfer function is given by [ a_0 + (a_1 * s) + (a_2 * s^2) +
/// (a_3 * s^3) + ... ] / [ b_0 + (b_1 * s) + (b_2 * s^2) + (b_3 * s^3) + ... ]
pub struct StateSpace {
    n_states: usize,
    state_values: Vec<f64>,
    a_vector: Vec<f64>,
    b_vector: Vec<f64>
}

impl StateSpace {
    pub fn new(a_vector: Vec<f64>, b_vector: Vec<f64>) -> Self {

        assert!(b_vector.len() >= a_vector.len());

        StateSpace {
            n_states: b_vector.len()-1,
            state_values: vec![0.0; b_vector.len()-1],
            a_vector: a_vector.into_iter().rev().collect(),
            b_vector: b_vector.into_iter().rev().collect()
        }
    }

    pub fn step(&mut self, timestep: f64, u: f64) -> f64 {
        let mut new_state_values = self.state_values.clone();

        let b_n = self.b_vector.iter().last().unwrap();
        let a_n = self.a_vector.get(self.b_vector.len()-1).unwrap_or(&0.0);
        
        let last_state_derivative: f64 = 
            ( u - (0..self.n_states).map(|i| self.state_values[i] * self.b_vector[i]).sum::<f64>() ) / self.b_vector.iter().last().unwrap();

        /*
        println!("Derivative: {}", last_state_derivative);
        println!("State: {:?}", self.state_values);
        */

        *new_state_values.iter_mut().last().unwrap() += last_state_derivative * timestep;

        for i in (0..(self.n_states-1)).rev() {
            new_state_values[i] += self.state_values[i+1] * timestep;
        }

        self.state_values = new_state_values;
        
        // self.state_values[0]

        let direct_response = (a_n / b_n) * u;

        let mut state_contribution = vec![0.0; self.b_vector.len()-1];
        for i in 0..self.b_vector.len()-1 {
            state_contribution[i] = (self.a_vector.get(i).unwrap_or(&0.0) - (a_n/b_n)*self.b_vector[i]) * self.state_values[i];
        }

        /*
        println!("State: {:?}", self.state_values);
        println!("Direct: {}", direct_response);
        println!("State contribution: {:?}", state_contribution);
        */

        return direct_response + state_contribution.iter().sum::<f64>()
    }
}


#[cfg(test)]
mod test {
    use super::*;

    mod state_space {
        use super::*;
        #[test]

        fn simple_response() {
            let mut state_space = StateSpace::new(vec![1.0], vec![1.0, 1.0]);

            for _ in 0..1000 {
                println!("{}", state_space.step(0.01, 1.0));
            }
        }
    }
}
