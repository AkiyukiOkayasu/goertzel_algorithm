//! Benchmark for the Goertzel algorithm.
//!
//! # Usage
//!
//! ```sh
//! cargo bench
//! ```
//!

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use goertzel_algorithm::{Goertzel, OptimizedGoertzel};

fn basic_goertzel_benchmark(c: &mut Criterion) {
    let mut goertzel = Goertzel::new();
    const SAMPLE_RATE: u32 = 48_000u32;
    const TARGET_FREQUENCY: f32 = 750.0f32;
    const BLOCK_SIZE: u32 = 128u32;
    goertzel.prepare(SAMPLE_RATE, TARGET_FREQUENCY, BLOCK_SIZE);
    c.bench_function("goertzel", |b| {
        b.iter(|| {
            for _i in 0..BLOCK_SIZE {
                let _ = goertzel.process_sample(black_box(&0.0f32));
            }
        })
    });
}

fn optimized_goertzel_benchmark(c: &mut Criterion) {
    let mut goertzel = OptimizedGoertzel::new();
    const SAMPLE_RATE: u32 = 48_000u32;
    const TARGET_FREQUENCY: f32 = 750.0f32;
    const BLOCK_SIZE: u32 = 128u32;
    goertzel.prepare(SAMPLE_RATE, TARGET_FREQUENCY, BLOCK_SIZE);
    c.bench_function("goertzel", |b| {
        b.iter(|| {
            for _i in 0..BLOCK_SIZE {
                let _ = goertzel.process_sample(black_box(&0.0f32));
            }
        })
    });
}

criterion_group!(
    benches,
    basic_goertzel_benchmark,
    optimized_goertzel_benchmark
);
criterion_main!(benches);
