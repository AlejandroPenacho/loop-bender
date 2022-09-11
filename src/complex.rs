#![allow(dead_code)]

use std::ops;

#[derive(Default,Clone,Copy)]
struct ComplexNumber {
    real: f64,
    imag: f64
}

impl ComplexNumber {
    fn from_cartesian(x: f64, y: f64) -> Self {
        ComplexNumber { real: x, imag: y }
    }

    fn from_mag_phase(mag: f64, phase: f64) -> Self {
        ComplexNumber { real: mag * phase.sin(), imag: mag * phase.cos() }
    }

    fn get_mag(&self) -> f64 {
        (self.real.powf(2.0) + self.imag.powf(2.0)).powf(0.5)
    }
    fn get_phase(&self) -> f64 {
        self.imag.atan2(self.real)
    }
}

impl ops::Add for ComplexNumber {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        ComplexNumber {
            real: self.real + other.real,
            imag: self.imag + other.imag
        }
    }
}

impl ops::Sub for ComplexNumber {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        ComplexNumber {
            real: self.real - other.real,
            imag: self.imag - other.imag
        }
    }

}

impl ops::Mul for ComplexNumber {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let mag = self.get_mag() * other.get_mag();
        let phase = self.get_phase() + other.get_phase();
        Self::from_mag_phase(mag, phase)
    }
}

impl ops::Div for ComplexNumber {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let mag = self.get_mag() / other.get_mag();
        let phase = self.get_phase() - other.get_phase();
        Self::from_mag_phase(mag, phase)
    }
}
