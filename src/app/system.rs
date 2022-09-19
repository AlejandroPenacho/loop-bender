use super::complex::ComplexNumber;

#[derive(Debug, Clone)]
pub enum PZElement {
    Pole(ComplexNumber),
    DoublePole(ComplexNumber),
    Zero(ComplexNumber),
    DoubleZero(ComplexNumber),
}

/// A system with dynamics
pub trait DynamicalSystem {
    fn get_steady_factor(&self) -> f64;
    fn get_pz_elements(&self) -> Vec<PZElement>;

    fn get_zeros(&self) -> Vec<ComplexNumber> {
        let mut output = Vec::new();
        for element in self.get_pz_elements() {
            match element {
                PZElement::Zero(z) => output.push(z),
                PZElement::DoubleZero(z) => {
                    output.push(z);
                    output.push(z.conjugate())
                }
                _ => {}
            }
        }
        output
    }

    fn get_poles(&self) -> Vec<ComplexNumber> {
        let mut output = Vec::new();
        for element in self.get_pz_elements() {
            match element {
                PZElement::Pole(z) => output.push(z),
                PZElement::DoublePole(z) => {
                    output.push(z);
                    output.push(z.conjugate())
                }
                _ => {}
            }
        }
        output
    }

    fn get_mag_at_freq(&self, freq: f64) -> f64 {
        let complex_freq = ComplexNumber::from_cartesian(0.0, freq);
        let mut magnitude = 1.0;

        for element in self.get_pz_elements() {
            match element {
                PZElement::Zero(zero) => {
                    magnitude *= (complex_freq - zero).get_mag();
                    if !zero.is_origin() {
                        magnitude /= zero.get_mag();
                    }
                }
                PZElement::Pole(pole) => {
                    magnitude /= (complex_freq - pole).get_mag();
                    if !pole.is_origin() {
                        magnitude *= pole.get_mag();
                    }
                }
                PZElement::DoubleZero(zero) => {
                    magnitude *= (complex_freq - zero).get_mag();
                    magnitude *= (complex_freq - zero.conjugate()).get_mag();

                    if !zero.is_origin() {
                        magnitude /= zero.get_mag().powf(2.0);
                    }
                }
                PZElement::DoublePole(pole) => {
                    magnitude /= (complex_freq - pole).get_mag();
                    magnitude /= (complex_freq - pole.conjugate()).get_mag();

                    if !pole.is_origin() {
                        magnitude *= pole.get_mag().powf(2.0);
                    }
                }
            }
        }

        magnitude * self.get_steady_factor()
    }

    fn get_phase_at_freq(&self, freq: f64) -> f64 {
        let complex_freq = ComplexNumber::from_cartesian(0.0, freq);
        let mut phase = 0.0;

        for element in self.get_pz_elements() {
            match element {
                PZElement::Zero(zero) => {
                    phase += (complex_freq - zero).get_phase();
                }
                PZElement::Pole(pole) => {
                    phase -= (complex_freq - pole).get_phase();
                }
                PZElement::DoubleZero(zero) => {
                    phase += (complex_freq - zero).get_phase();
                    phase += (complex_freq - zero.conjugate()).get_phase();
                }
                PZElement::DoublePole(pole) => {
                    phase -= (complex_freq - pole).get_phase();
                    phase -= (complex_freq - pole.conjugate()).get_phase();
                }
            }
        }
        phase
    }

    fn get_state_space(&self) -> StateSpace {
        let pz_elements = self.get_pz_elements();
        let mut a_vector = vec![1.0];
        let mut b_vector = vec![1.0];

        for pz in pz_elements {
            match pz {
                PZElement::Zero(z) => {
                    a_vector = polynomial_product(&a_vector, &[-z.as_cartesian()[0], 1.0]);
                }
                PZElement::Pole(p) => {
                    b_vector = polynomial_product(&b_vector, &[-p.as_cartesian()[0], 1.0]);
                }
                PZElement::DoubleZero(zz) => {
                    let real = zz.as_cartesian()[0];
                    let imag = zz.as_cartesian()[1];
                    a_vector = polynomial_product(
                        &a_vector,
                        &[real.powf(2.0) + imag.powf(2.0), -2.0 * real, 1.0],
                    );
                }
                PZElement::DoublePole(pp) => {
                    let real = pp.as_cartesian()[0];
                    let imag = pp.as_cartesian()[1];
                    b_vector = polynomial_product(
                        &b_vector,
                        &[real.powf(2.0) + imag.powf(2.0), -2.0 * real, 1.0],
                    );
                }
            }
        }

        let polynomical_gain = a_vector.iter().find(|&&x| x > 0.00001).unwrap()
            / b_vector.iter().find(|&&x| x > 0.00001).unwrap();

        a_vector
            .iter_mut()
            .for_each(|x| *x *= self.get_steady_factor() / polynomical_gain);

        StateSpace::new(a_vector, b_vector)
    }
}

fn polynomial_product(p_1: &[f64], p_2: &[f64]) -> Vec<f64> {
    let mut product = vec![0.0; p_1.len() + p_2.len() - 1];
    for (i, value_a) in p_1.iter().enumerate() {
        for (j, value_b) in p_2.iter().enumerate() {
            product[i + j] += value_a * value_b;
        }
    }
    product
}

impl<T: DynamicalSystem> DynamicalSystem for &[T] {
    fn get_steady_factor(&self) -> f64 {
        self.iter().map(|x| x.get_steady_factor()).product::<f64>()
    }
    fn get_pz_elements(&self) -> Vec<PZElement> {
        self.iter().flat_map(|x| x.get_pz_elements()).collect()
    }
}

#[derive(Default)]
pub struct Controller {
    pub pid: Pid,
    pub lead_lag: LLController,
    pub lag_lead: LLController,
}

impl DynamicalSystem for Controller {
    fn get_steady_factor(&self) -> f64 {
        self.pid.get_steady_factor()
            * self.lead_lag.get_steady_factor()
            * self.lag_lead.get_steady_factor()
    }

    fn get_pz_elements(&self) -> Vec<PZElement> {
        let mut output = self.pid.get_pz_elements();
        output.append(&mut self.lead_lag.get_pz_elements());
        output.append(&mut self.lag_lead.get_pz_elements());

        output
    }
}

/// A Pid controller
pub struct Pid {
    k: f64,
    t_i: f64,
    t_d: f64,
    tau: f64,
    integral_on: bool,
    derivative_on: bool,
}

impl Pid {
    pub fn get_k(&mut self) -> &mut f64 {
        &mut self.k
    }
    pub fn get_ti(&mut self) -> &mut f64 {
        &mut self.t_i
    }
    pub fn get_td(&mut self) -> &mut f64 {
        &mut self.t_d
    }
    pub fn get_integral_on(&mut self) -> &mut bool {
        &mut self.integral_on
    }
    pub fn get_derivative_on(&mut self) -> &mut bool {
        &mut self.derivative_on
    }
}

impl DynamicalSystem for Pid {
    fn get_steady_factor(&self) -> f64 {
        self.k
    }

    /*
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
    */

    fn get_pz_elements(&self) -> Vec<PZElement> {
        let mut output = Vec::new();
        if self.integral_on {
            output.push(PZElement::Pole(ComplexNumber::from_cartesian(0.0, 0.0)));
        }
        if self.derivative_on {
            output.push(PZElement::Pole(ComplexNumber::from_cartesian(
                -1.0 / self.tau,
                0.0,
            )));
        }

        if self.integral_on && !self.derivative_on {
            output.push(PZElement::Zero(ComplexNumber::from_cartesian(
                -1.0 / self.t_i,
                0.0,
            )));
        } else if !self.integral_on && self.derivative_on {
            output.push(PZElement::Zero(ComplexNumber::from_cartesian(
                -1.0 / (self.t_d + self.tau),
                0.0,
            )));
        } else if self.integral_on && self.derivative_on {
            let midpoint =
                -((self.t_i + self.tau) / (2.0 * self.t_i * (self.t_d + self.tau).powf(2.0)));

            let delta = ((self.t_i + self.tau).powf(2.0)
                / (4.0 * self.t_i.powf(2.0) * (self.t_d + self.tau).powf(2.0)))
                - 1.0 / (self.t_i * (self.t_d + self.tau));

            if delta >= 0.0 {
                output.push(PZElement::Zero(ComplexNumber::from_cartesian(
                    midpoint + delta.powf(0.5),
                    0.0,
                )));
                output.push(PZElement::Zero(ComplexNumber::from_cartesian(
                    midpoint - delta.powf(0.5),
                    0.0,
                )));
            } else {
                output.push(PZElement::DoubleZero(ComplexNumber::from_cartesian(
                    midpoint,
                    (-delta).powf(0.5),
                )));
            }
        }

        output
    }
}

impl Default for Pid {
    fn default() -> Self {
        Pid {
            k: 1.0,
            t_i: 20.0,
            t_d: 1.0,
            tau: 0.05,
            integral_on: false,
            derivative_on: false,
        }
    }
}

pub struct LLController {
    zero_frequency: f64,
    pole_frequency: f64,
    activated: bool,
}

impl Default for LLController {
    fn default() -> Self {
        Self {
            zero_frequency: 1.0,
            pole_frequency: 1.0,
            activated: false,
        }
    }
}

impl LLController {
    pub fn is_activated(&mut self) -> &mut bool {
        &mut self.activated
    }
    pub fn get_zero_freq(&mut self) -> &mut f64 {
        &mut self.zero_frequency
    }
    pub fn get_pole_freq(&mut self) -> &mut f64 {
        &mut self.pole_frequency
    }
}

impl DynamicalSystem for LLController {
    fn get_steady_factor(&self) -> f64 {
        1.0
    }

    fn get_pz_elements(&self) -> Vec<PZElement> {
        if self.activated {
            vec![
                PZElement::Zero(ComplexNumber::from_cartesian(-self.zero_frequency, 0.0)),
                PZElement::Pole(ComplexNumber::from_cartesian(-self.pole_frequency, 0.0)),
            ]
        } else {
            vec![]
        }
    }
}

/// A generic dynamic system, defined by its gain, poles and zeros
#[derive(Clone)]
pub struct Model {
    components: Vec<PZElement>,
    gain: f64,
}

impl Default for Model {
    fn default() -> Model {
        let components = vec![
            PZElement::Pole(ComplexNumber::from_cartesian(-2.0, 0.0)),
            PZElement::Pole(ComplexNumber::from_cartesian(-3.0, 0.0)),
        ];
        Model {
            components,
            gain: 1.0,
        }
    }
}

impl Model {
    pub fn push_pole(&mut self, coord: (f64, f64)) {
        if coord.1.abs() < coord.0.abs() / 10.0 {
            self.components
                .push(PZElement::Pole(ComplexNumber::from_cartesian(coord.0, 0.0)));
        } else {
            self.components
                .push(PZElement::DoublePole(ComplexNumber::from_cartesian(
                    coord.0, coord.1,
                )));
        }
    }

    pub fn push_zero(&mut self, coord: (f64, f64)) {
        if coord.1.abs() < coord.0.abs() / 10.0 {
            self.components
                .push(PZElement::Zero(ComplexNumber::from_cartesian(coord.0, 0.0)));
        } else {
            self.components
                .push(PZElement::DoubleZero(ComplexNumber::from_cartesian(
                    coord.0, coord.1,
                )));
        }
    }

    pub fn remove_closest_element(&mut self, coord: (f64, f64)) {
        let x = coord.0;
        let y = coord.1;

        let closest_element = self
            .components
            .iter()
            .map(|element| match element {
                PZElement::Zero(point) | PZElement::Pole(point) => {
                    let [p_x, p_y] = point.as_cartesian();
                    (p_x - x).powf(2.0) + (p_y - y).powf(2.0)
                }
                PZElement::DoubleZero(point) | PZElement::DoublePole(point) => {
                    let [p_x, p_y] = point.as_cartesian();
                    (p_x - x).powf(2.0) + (p_y - y).powf(2.0).min((p_y + y).powf(2.0))
                }
            })
            .enumerate()
            .reduce(|(i_min, dist_min), (i_new, dist_new)| {
                if let Some(std::cmp::Ordering::Less) = dist_new.partial_cmp(&dist_min) {
                    (i_new, dist_new)
                } else {
                    (i_min, dist_min)
                }
            })
            .map(|x| x.0);

        if let Some(index) = closest_element {
            self.components.remove(index);
        }
    }

    pub fn link_system<T: DynamicalSystem>(&self, rhs: &T) -> Model {
        let mut all_elements = self.get_pz_elements();
        let mut gain = self.get_steady_factor();

        all_elements.append(&mut rhs.get_pz_elements());
        gain *= rhs.get_steady_factor();

        Model {
            gain,
            components: all_elements,
        }
    }
}

impl DynamicalSystem for Model {
    fn get_steady_factor(&self) -> f64 {
        self.gain
    }

    fn get_pz_elements(&self) -> Vec<PZElement> {
        self.components.clone()
    }
}

/// Represents a set of differential equations that allow the simulation of a
/// dynamic system
///
/// The transfer function is given by [ a_0 + (a_1 * s) + (a_2 * s^2) +
/// (a_3 * s^3) + ... ] / [ b_0 + (b_1 * s) + (b_2 * s^2) + (b_3 * s^3) + ... ]
pub struct StateSpace {
    n_states: usize,
    last_input: f64,
    state_values: Vec<f64>,
    pub a_vector: Vec<f64>,
    pub b_vector: Vec<f64>,
}

impl StateSpace {
    pub fn new(a_vector: Vec<f64>, b_vector: Vec<f64>) -> Self {
        assert!(b_vector.len() >= a_vector.len());

        StateSpace {
            n_states: b_vector.len() - 1,
            last_input: 0.0,
            state_values: vec![0.0; b_vector.len() - 1],
            a_vector,
            b_vector,
        }
    }

    pub fn get_output(&self) -> f64 {
        let b_n = self.b_vector.iter().last().unwrap();
        let a_n = self.a_vector.get(self.b_vector.len() - 1).unwrap_or(&0.0);

        if self.n_states == 0 {
            return self.last_input * a_n / b_n;
        }

        let direct_response = (a_n / b_n) * self.last_input;

        let mut state_contribution = vec![0.0; self.b_vector.len() - 1];
        for (i, b_i) in self
            .b_vector
            .iter()
            .enumerate()
            .take(self.b_vector.len() - 1)
        {
            state_contribution[i] =
                (self.a_vector.get(i).unwrap_or(&0.0) - (a_n / b_n) * b_i) * self.state_values[i];
        }

        direct_response + state_contribution.iter().sum::<f64>()
    }

    pub fn step(&mut self, timestep: f64, u: f64) {
        self.last_input = u;
        if self.n_states == 0 {
            return;
        }

        let mut new_state_values = self.state_values.clone();

        /*
        let b_n = self.b_vector.iter().last().unwrap();
        let a_n = self.a_vector.get(self.b_vector.len() - 1).unwrap_or(&0.0);
        */

        let last_state_derivative: f64 = (u
            - (0..self.n_states)
                .map(|i| self.state_values[i] * self.b_vector[i])
                .sum::<f64>())
            / self.b_vector.iter().last().unwrap();

        /*
        println!("Derivative: {}", last_state_derivative);
        println!("State: {:?}", self.state_values);
        */

        if let Some(value) = new_state_values.iter_mut().last() {
            *value += last_state_derivative * timestep
        }

        for i in (0..(self.n_states - 1)).rev() {
            new_state_values[i] += self.state_values[i + 1] * timestep;
        }

        self.state_values = new_state_values;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn polynomal_multiplication_1() {
        //
        //  (2 + 3*s + 5*s^2) * (-3 + 2*s) =
        //  (-6 - 9*s - 15*s^2) + (4*s + 6*s^2+ 10*s^3) =
        //  (-6 - 5*s - 9*s^2 + 10*s^3)
        //

        let expected = [-6.0, -5.0, -9.0, 10.0];
        let product = polynomial_product(&[2.0, 3.0, 5.0], &[-3.0, 2.0]);

        for (x, y) in expected.iter().zip(product.iter()) {
            assert!((x - y).abs() < 0.00001);
        }
    }

    mod state_space {
        use super::*;
        #[test]

        fn simple_response() {
            let mut state_space = StateSpace::new(vec![1.0], vec![1.0, 1.0]);

            for _ in 0..1000 {
                state_space.step(0.01, 1.0);
            }
        }
    }
}
