# Goertzel algorithm  

[![Cargo](https://img.shields.io/crates/v/goertzel_algorithm.svg)](https://crates.io/crates/goertzel_algorithm)
[![Documentation](https://docs.rs/goertzel_algorithm/badge.svg)](https://docs.rs/goertzel_algorithm)

Useful when analyzing the amplitude or phase of a specific frequency.  

## Difference from FFT

When analyzing only a few specific frequencies, it may be more efficient than an FFT.  
Different from the FFT, the computational cost remains the same even when the block size is not a power of 2.  

## no_std

Works with no_std by default.
