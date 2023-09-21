//! # Goertzel algorithm
//!
//! This crate provides a Goertzel algorithm implementation.
//!
//! # Example
//!
//! ```
//! use goertzel_algorithm::Goertzel;
//! use approx::{assert_relative_eq, assert_relative_ne};
//!
//! const SAMPLE_RATE: u32 = 48_000u32;
//! const TARGET_FREQUENCY: f32 = 750.0f32;
//! const BLOCK_SIZE: u32 = 128u32;
//! let phase_increment = TARGET_FREQUENCY * std::f32::consts::PI * 2.0f32 * (1.0f32 / SAMPLE_RATE as f32);
//! let mut phase = 0.0f32;
//!
//! let mut goertzel = Goertzel::new();
//! goertzel.prepare(SAMPLE_RATE, TARGET_FREQUENCY, BLOCK_SIZE);
//!
//! for i in 0..BLOCK_SIZE {
//!     let input = phase.sin();//Generate a sine wave same frequency as the target frequency
//!     if let Some(mag_phase) = goertzel.process_sample(&input) {
//!         println!("{}: {}", i, mag_phase.magnitude);//127: 1.0
//!     }
//!
//!     phase += phase_increment;
//!     if phase >= std::f32::consts::PI * 2.0f32 {
//!         phase -= std::f32::consts::PI * 2.0f32;
//!     }
//! }            
//! ```

#![cfg_attr(not(test), no_std)]

use core::f32::consts::PI;

#[derive(Default)]
pub struct MagnitudePhase {
    pub magnitude: f32,
    pub phase: f32,
}

/// Basic Goertzel algorithm.
///
/// It calculates the magnitude and phase of the frequency.
#[derive(Default)]
pub struct Goertzel {
    sample_rate: u32,
    target_frequency: f32,
    q0: f32,
    q1: f32,
    q2: f32,
    coeff: f32,
    cosine: f32,
    sine: f32,
    block_size: u32,
    counter: u32,
}

impl Goertzel {
    pub fn new() -> Self {
        Self {
            sample_rate: 48_000u32,
            target_frequency: 0.0f32,
            q0: 0.0f32,
            q1: 0.0f32,
            q2: 0.0f32,
            coeff: 0.0f32,
            cosine: 0.0f32,
            sine: 0.0f32,
            block_size: 128u32,
            counter: 0u32,
        }
    }

    /// Prepare the Goertzel algorithm's coefficients.
    ///
    /// # Arguments
    ///
    /// * 'sample_rate' - The sample rate of the input signal.
    /// * 'target_frequency' - The frequency to detect.
    /// * 'block_size' - The number of samples to process at a time. It is like FFT size, but does not have to be a power of 2.
    pub fn prepare(&mut self, sample_rate: u32, target_frequency: f32, block_size: u32) {
        let k = (block_size as f32 * target_frequency) / sample_rate as f32;
        let w = (2.0f32 * PI / block_size as f32) * k;
        self.cosine = libm::cosf(w);
        self.sine = libm::sinf(w);
        self.coeff = 2.0f32 * self.cosine;

        self.sample_rate = sample_rate;
        self.target_frequency = target_frequency;
        self.block_size = block_size;

        self.q0 = 0.0f32;
        self.q1 = 0.0f32;
        self.q2 = 0.0f32;
        self.counter = 0u32;
    }

    /// Process a sample.
    ///
    /// # Arguments
    ///
    /// * 'input' - The sample to add.
    ///
    /// # Returns
    ///
    /// Returns the magnitude and phase of the frequency if a block is completed, otherwise None.
    /// Magnitude is normalized to [0.0, 1.0].
    #[inline]
    #[must_use]
    pub fn process_sample(&mut self, input: &f32) -> Option<MagnitudePhase> {
        self.q0 = self.coeff * self.q1 - self.q2 + input;
        self.q2 = self.q1;
        self.q1 = self.q0;

        self.counter += 1;
        if self.counter >= self.block_size {
            // Basic goertzel (magnitude and phase)
            let real = self.q1 - self.q2 * self.cosine;
            let imag = self.q2 * self.sine;
            let magnitude = libm::sqrtf(real * real + imag * imag);
            let phase = libm::atan2f(imag, real);

            // optimized goertzel (only magnitude)
            /*
            let magnitude = libm::sqrtf(
                (self.q1 * self.q1) + (self.q2 * self.q2) - (self.q1 * self.q2 * self.coeff),
            );
            */

            let normalized_magnitude = magnitude / (self.block_size as f32 / 2.0f32); //[0.0, 1.0]

            self.q1 = 0.0f32;
            self.q2 = 0.0f32;
            self.counter = 0u32;

            Some(MagnitudePhase {
                magnitude: normalized_magnitude,
                phase,
            })
        } else {
            None
        }
    }
}

/// Optimized Goertzel algorithm.
///
/// It only calculates the magnitude.
/// ```
/// use goertzel_algorithm::OptimizedGoertzel;
/// use approx::{assert_relative_eq, assert_relative_ne};
///
/// const SAMPLE_RATE: u32 = 48_000u32;
/// const TARGET_FREQUENCY: f32 = 750.0f32;
/// const BLOCK_SIZE: u32 = 128u32;
/// let phase_increment = TARGET_FREQUENCY * std::f32::consts::PI * 2.0f32 * (1.0f32 / SAMPLE_RATE as f32);
/// let mut phase = 0.0f32;
///
/// let mut goertzel = OptimizedGoertzel::new();
/// goertzel.prepare(SAMPLE_RATE, TARGET_FREQUENCY, BLOCK_SIZE);
///
/// for i in 0..BLOCK_SIZE {
///     let input = phase.sin();//Generate a sine wave same frequency as the target frequency
///     if let Some(magnitude) = goertzel.process_sample(&input) {
///         println!("{}: {}", i, magnitude);//127: 1.0
///     }
///
///     phase += phase_increment;
///     if phase >= std::f32::consts::PI * 2.0f32 {
///         phase -= std::f32::consts::PI * 2.0f32;
///     }
/// }            
/// ```
#[derive(Default)]
pub struct OptimizedGoertzel {
    sample_rate: u32,
    target_frequency: f32,
    q0: f32,
    q1: f32,
    q2: f32,
    coeff: f32,
    block_size: u32,
    counter: u32,
}

impl OptimizedGoertzel {
    pub fn new() -> Self {
        Self {
            sample_rate: 48_000u32,
            target_frequency: 0.0f32,
            q0: 0.0f32,
            q1: 0.0f32,
            q2: 0.0f32,
            coeff: 0.0f32,
            block_size: 128u32,
            counter: 0u32,
        }
    }

    /// Prepare the Goertzel algorithm's coefficients.
    ///
    /// # Arguments
    ///
    /// * 'sample_rate' - The sample rate of the input signal.
    /// * 'target_frequency' - The frequency to detect.
    /// * 'block_size' - The number of samples to process at a time. It is like FFT size, but does not have to be a power of 2.
    pub fn prepare(&mut self, sample_rate: u32, target_frequency: f32, block_size: u32) {
        let k = (block_size as f32 * target_frequency) / sample_rate as f32;
        let w = (2.0f32 * PI / block_size as f32) * k;
        self.coeff = 2.0f32 * libm::cosf(w);

        self.sample_rate = sample_rate;
        self.target_frequency = target_frequency;
        self.block_size = block_size;

        self.q0 = 0.0f32;
        self.q1 = 0.0f32;
        self.q2 = 0.0f32;
        self.counter = 0u32;
    }

    /// Process a sample.
    ///
    /// # Arguments
    ///
    /// * 'input' - The sample to add.
    ///
    /// # Returns
    ///
    /// Returns the magnitude of the frequency if a block is completed, otherwise None.
    /// Magnitude is normalized to [0.0, 1.0].
    #[inline]
    #[must_use]
    pub fn process_sample(&mut self, input: &f32) -> Option<f32> {
        self.q0 = self.coeff * self.q1 - self.q2 + input;
        self.q2 = self.q1;
        self.q1 = self.q0;

        self.counter += 1;
        if self.counter >= self.block_size {
            // optimized goertzel (only magnitude)
            let magnitude = libm::sqrtf(
                (self.q1 * self.q1) + (self.q2 * self.q2) - (self.q1 * self.q2 * self.coeff),
            );

            let normalized_magnitude = magnitude / (self.block_size as f32 / 2.0f32); //[0.0, 1.0]

            self.q1 = 0.0f32;
            self.q2 = 0.0f32;
            self.counter = 0u32;

            Some(normalized_magnitude)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_relative_eq, assert_relative_ne};

    #[test]
    fn same_freq_basic() {
        const SAMPLE_RATE: u32 = 48_000u32;
        const TARGET_FREQUENCY: f32 = 750.0f32;
        const BLOCK_SIZE: u32 = 128u32;
        let phase_increment =
            TARGET_FREQUENCY * std::f32::consts::PI * 2.0f32 * (1.0f32 / SAMPLE_RATE as f32);
        let mut phase = 0.0f32;

        let mut goertzel = Goertzel::new();
        goertzel.prepare(SAMPLE_RATE, TARGET_FREQUENCY, BLOCK_SIZE);
        for _k in 0..100 {
            for i in 0..BLOCK_SIZE {
                let input = phase.sin();
                if let Some(v) = goertzel.process_sample(&input) {
                    assert_eq!(i, BLOCK_SIZE - 1);
                    assert_relative_eq!(v.magnitude, 1.0f32, epsilon = 0.0001f32);
                }
                phase += phase_increment;
                if phase >= std::f32::consts::PI * 2.0f32 {
                    phase -= std::f32::consts::PI * 2.0f32;
                }
            }
        }
    }

    #[test]
    fn same_freq_optimized() {
        const SAMPLE_RATE: u32 = 48_000u32;
        const TARGET_FREQUENCY: f32 = 750.0f32;
        const BLOCK_SIZE: u32 = 128u32;
        let phase_increment =
            TARGET_FREQUENCY * std::f32::consts::PI * 2.0f32 * (1.0f32 / SAMPLE_RATE as f32);
        let mut phase = 0.0f32;

        let mut goertzel = OptimizedGoertzel::new();
        goertzel.prepare(SAMPLE_RATE, TARGET_FREQUENCY, BLOCK_SIZE);
        for _k in 0..100 {
            for i in 0..BLOCK_SIZE {
                let input = phase.sin();
                if let Some(mag) = goertzel.process_sample(&input) {
                    assert_eq!(i, BLOCK_SIZE - 1);
                    assert_relative_eq!(mag, 1.0f32, epsilon = 0.0001f32);
                }
                phase += phase_increment;
                if phase >= std::f32::consts::PI * 2.0f32 {
                    phase -= std::f32::consts::PI * 2.0f32;
                }
            }
        }
    }

    #[test]
    fn different_freq_basic() {
        const SAMPLE_RATE: u32 = 48_000u32;
        const TARGET_FREQUENCY: f32 = 750.0f32;
        const BLOCK_SIZE: u32 = 128u32;
        let phase_increment =
            TARGET_FREQUENCY * std::f32::consts::PI * 2.0f32 * (1.0f32 / SAMPLE_RATE as f32);
        let mut phase = 0.0f32;

        let mut goertzel = Goertzel::new();
        goertzel.prepare(SAMPLE_RATE, 2000.0f32, BLOCK_SIZE);

        for _k in 0..100 {
            for i in 0..BLOCK_SIZE {
                let input = phase.sin();
                if let Some(v) = goertzel.process_sample(&input) {
                    assert_eq!(i, BLOCK_SIZE - 1);
                    assert_relative_ne!(v.magnitude, 0.0001f32);
                }
                phase += phase_increment;
                if phase >= std::f32::consts::PI * 2.0f32 {
                    phase -= std::f32::consts::PI * 2.0f32;
                }
            }
        }
    }

    #[test]
    fn different_freq_optimized() {
        const SAMPLE_RATE: u32 = 48_000u32;
        const TARGET_FREQUENCY: f32 = 750.0f32;
        const BLOCK_SIZE: u32 = 128u32;
        let phase_increment =
            TARGET_FREQUENCY * std::f32::consts::PI * 2.0f32 * (1.0f32 / SAMPLE_RATE as f32);
        let mut phase = 0.0f32;

        let mut goertzel = OptimizedGoertzel::new();
        goertzel.prepare(SAMPLE_RATE, 2000.0f32, BLOCK_SIZE);

        for _k in 0..100 {
            for i in 0..BLOCK_SIZE {
                let input = phase.sin();
                if let Some(mag) = goertzel.process_sample(&input) {
                    assert_eq!(i, BLOCK_SIZE - 1);
                    assert_relative_ne!(mag, 0.0001f32);
                }
                phase += phase_increment;
                if phase >= std::f32::consts::PI * 2.0f32 {
                    phase -= std::f32::consts::PI * 2.0f32;
                }
            }
        }
    }
}
