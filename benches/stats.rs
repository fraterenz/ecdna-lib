use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ecdna_lib::distribution::EcDNADistribution;
use rand::Rng;

fn compute_variance(c: &mut Criterion) {
    static KB: usize = 1024;

    let mut group = c.benchmark_group("compute_variance");
    let mut rng = rand::thread_rng();

    for size in [KB, 2 * KB, 4 * KB, 8 * KB, 16 * KB].iter() {
        group.throughput(Throughput::Bytes(*size as u64));

        let my_vec: Vec<u16> = (0..*size).map(|_| rng.gen_range(0..u16::MAX)).collect();
        let distribution = EcDNADistribution::from(my_vec);

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| {
                distribution.compute_variance();
            })
        });
    }

    group.finish();
}

fn compute_mean(c: &mut Criterion) {
    static KB: usize = 1024;

    let mut group = c.benchmark_group("compute_mean");
    let mut rng = rand::thread_rng();

    for size in [KB, 2 * KB, 4 * KB, 8 * KB, 16 * KB].iter() {
        group.throughput(Throughput::Bytes(*size as u64));
        let my_vec: Vec<u16> = (0..*size).map(|_| rng.gen_range(0..u16::MAX)).collect();
        let distribution = EcDNADistribution::from(my_vec);

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| {
                distribution.compute_mean();
            })
        });
    }

    group.finish();
}

criterion_group!(benches, compute_mean, compute_variance);
criterion_main!(benches);
