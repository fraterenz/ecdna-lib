# Changelog

## v0.1.1
- Use the Kolmogorov-Smirnov implementation found in [`kolmogorov_smirnov`](https://github.com/daithiocrualaoich/kolmogorov_smirnov/blob/cb067e92ec837efbad66e8bbcf85500ad778feb8/src/test.rs#L127) crate.
- Run the ks test considering only cells with ecDNAs.

## v0.1.2
- New method `EcDNADistribution::drop_cells_with_k_copies`
