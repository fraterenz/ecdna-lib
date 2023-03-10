# Changelog

## v0.1.1
- Use the Kolmogorov-Smirnov implementation found in [`kolmogorov_smirnov`](https://github.com/daithiocrualaoich/kolmogorov_smirnov/blob/cb067e92ec837efbad66e8bbcf85500ad778feb8/src/test.rs#L127) crate.
- Run the ks test considering only cells with ecDNAs.

## v0.1.2
- New method `EcDNADistribution::drop_cells_with_k_copies`

## v0.1.3
### Bugfix
Fix the method `EcDNADistribution::drop_cells_with_k_copies`.

## v0.2.0
Remove the rng `ChaCha` and keep it generic with `impl Rng` trait.

## v0.2.2
Make `distribution.create_histogram` public and remove `create_histogram_f32`.

## v0.3.0
Make `undersample` private and `sample` public which samples the ecDNA distribution according to a `SamplingStrategy`.

## v0.3.1
Allows `sample` to take a sample with a number of cells equal to the ecDNA distribution size, that is to not undersample.

## v0.3.2
Implement `Eq` for `SamplingStrategy`.

## v0.3.4
Implement `SamplingStrategy::Poisson`.
### Bugfix
Ensure that cells without ecDNAs cannot gain any ecDNA by randomly mapping the distribution to a Gaussian or Poisson.

## v0.3.5
Implement `SamplingStrategy::Exponential(scale)`. Note that this distribution has the same scale parameters for all cells, that is independently of the number of `k` copies present in cells.

