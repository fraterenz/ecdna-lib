use crate::DNACopy;
use anyhow::{bail, Context};
use rand::{seq::SliceRandom, Rng};
use rand_distr::Distribution;
use std::cmp::{min, Ord};
use std::fmt;
use std::{collections::HashMap, fs, num::NonZeroU16, path::Path};

/// Sampling strategies to sample the [`EcDNADistribution`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SamplingStrategy {
    /// Randomly pick cells from the ecDNA distribution.
    Uniform,
    /// Map each entry of the ecDNA distribution `k` into a Normal distribution
    /// with mean `k` and std of 1, and then randomly pick cells from this
    /// modified distribution.
    ///
    /// Note that if the distribution contains only cells w/o any ecDNAs, then
    /// this will not generate any new cell with ecDNA (that is it's a biased
    /// Gaussian).
    Gaussian(Sigma),
    /// Map each entry of the ecDNA distribution `k` into a Poisson point
    /// process with mean `lambda * k` and , and then randomly pick cells from
    /// this modified distribution.
    Poisson(Lambda),
    /// Map each entry of the ecDNA distribution `k` into a Exponential
    /// distribution with a rate parameter, and then randomly pick cells from
    /// this modified distribution.
    ///
    /// The mapping is `k * (1-exp(rate))`, that is then ones we have before
    /// transcription `k` times the ones that have not degraded `1-exp(rate)`.
    Exponential(Lambda),
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// The standard deviation of the Normal distribution used to sample the ecDNA
/// distribution with [`SamplingStrategy::Gaussian`].
///
/// Is a float that cannot be smaller or equal than 0 nor infinite.
pub struct Sigma(Lambda);

impl Sigma {
    pub fn new(sigma: f32) -> Sigma {
        //! ## Panics
        //! When sigma is smaller or equal than 0 or is not finite.
        Sigma(Lambda::new(sigma))
    }

    pub fn get(&self) -> f32 {
        self.0 .0
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// The lambda of the Poisson point process used to sample the ecDNA
/// distribution with [`SamplingStrategy::Poisson`].
///
/// Is a float that cannot be smaller or equal than 0 nor infinite.
pub struct Lambda(f32);

impl Lambda {
    pub fn new(lambda: f32) -> Lambda {
        //! ## Panics
        //! When lambda is smaller or equal than 0 or is not finite.
        assert!(
            lambda.is_normal() && lambda.is_sign_positive(),
            "lambda {} not valid!, must be finite and > 0",
            lambda
        );
        Lambda(lambda)
    }

    pub fn get(&self) -> f32 {
        self.0
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct EcDNADistribution {
    nminus: u64,
    nplus: Vec<DNACopy>,
}

impl fmt::Display for EcDNADistribution {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut hist = self
            .create_histogram()
            .into_iter()
            .collect::<Vec<(u16, u64)>>();
        hist.sort_by_key(|(ecdna, _)| *ecdna);
        hist.into_iter()
            .try_for_each(|(ecdna, cells)| writeln!(f, "{}:{}", ecdna, cells))
    }
}

impl From<Vec<u16>> for EcDNADistribution {
    fn from(vec: Vec<u16>) -> Self {
        let mut nplus = vec;
        nplus.shrink_to_fit();
        nplus.sort();
        let (nminus, nplus) = nplus.split_at(nplus.partition_point(|ele| ele == &0u16));
        Self {
            nminus: nminus.len() as u64,
            nplus: nplus
                .iter()
                .map(|&k| unsafe { NonZeroU16::new_unchecked(k) })
                .collect(),
        }
    }
}

impl EcDNADistribution {
    pub fn new(distribution: HashMap<u16, u64>, iterations: usize) -> Self {
        let nminus = *distribution.get(&0).unwrap_or(&0);

        let mut nplus = Vec::with_capacity(iterations);
        for (copies, cells) in distribution.into_iter() {
            if copies > 0 {
                for _ in 0..cells {
                    nplus.push(unsafe { NonZeroU16::new_unchecked(copies) });
                }
            }
        }

        EcDNADistribution { nminus, nplus }
    }

    pub fn get_nminus(&self) -> &u64 {
        &self.nminus
    }

    pub fn is_nplus_empty(&self) -> bool {
        //! No cells with ecDNAs are left in the population.
        self.nplus.is_empty()
    }

    pub fn pick_remove_random_nplus(&mut self, rng: &mut impl Rng) -> anyhow::Result<DNACopy> {
        //! Returns ecDNA copies of a nplus cell and remove it from the ecDNA
        //! distribution.
        //! Note that all cells with ecDNAs have the same probability of being
        //! choosen, which is true only when we consider constant fitness.
        //!
        //! ## Fails
        //! Fails when there are no cells with ecDNA in the population, that is
        //! [`EcDNADistribution::is_nplus_empty`].
        if self.is_nplus_empty() {
            bail!("Empty distr")
        } else {
            let idx = self.pick_random_nplus(rng);
            Ok(self.nplus.swap_remove(idx))
        }
    }

    fn pick_random_nplus(&self, rng: &mut impl Rng) -> usize {
        rng.gen_range(0..self.nplus.len())
    }

    pub fn increase_nplus(&mut self, copies: Vec<DNACopy>, verbosity: u8) {
        //! Update the ecdna distribution by adding `copies` to the population
        //! of cells with ecDNAs.
        for k in copies {
            if verbosity > 1 {
                println!("Increasing count for clone {}", k);
            }
            self.nplus.push(k);
        }
    }

    pub fn increase_nminus(&mut self) {
        self.nminus += 1;
    }

    pub fn decrease_nplus(&mut self, rng: &mut impl Rng, verbosity: u8) {
        let idx = self.pick_random_nplus(rng);

        if verbosity > 1 {
            println!("Clone {} will loose one cell", self.nplus[idx]);
        }
        self.nplus.swap_remove(idx);
    }

    pub fn decrease_nminus(&mut self) {
        self.nminus -= 1;
    }

    pub fn compute_nplus(&self) -> u64 {
        self.nplus.len() as u64
    }

    pub fn is_empty(&self) -> bool {
        self.is_nplus_empty() && self.nminus == 0
    }

    pub fn load(path2file: &Path, capacity: usize) -> anyhow::Result<Self> {
        let extension = path2file
            .extension()
            .with_context(|| format!("Do not recognize extension of file {:#?}", path2file))
            .unwrap();
        if let Some("json") = extension.to_str() {
            let path2read = fs::read_to_string(path2file)
                .with_context(|| format!("Cannot read {:#?}", path2file))?;

            let map: HashMap<u16, u64> = serde_json::from_str(&path2read)
                .map_err(|e| anyhow::anyhow!(e))
                .with_context(|| "Cannot load ecDNA distribution")?;

            return Ok(EcDNADistribution::new(map, capacity));
        }
        bail!("Extension not recognized, must be JSON `.json`")
    }

    pub fn save(&self, path2file: &Path, verbosity: u8) -> anyhow::Result<()> {
        let map: HashMap<u16, u64> = self.create_histogram();
        let distribution =
            serde_json::to_string(&map).expect("Cannot serialize the ecDNA distribution");
        if verbosity > 1 {
            println!(
                "Saving ecDNA distribution to {:#?}\n{:#?}",
                path2file, distribution
            );
        }
        fs::create_dir_all(path2file.parent().unwrap()).expect("Cannot create dir");
        fs::write(path2file, distribution)
            .with_context(|| "Cannot save the ecDNA distribution".to_string())?;
        Ok(())
    }

    pub fn sample<R: Rng + Clone>(&mut self, nb_cells: u64, strategy: &SamplingStrategy, rng: R) {
        //! Draw a random sample without replacement from the
        //! `EcDNADistribution` according to the [`SamplingStrategy`].
        //!
        //! ## Panics
        //! Panics when `nb_cells` is greater than the cells in the
        //! distribution or when `nb_cells` is 0.
        assert!(nb_cells <= (*self.get_nminus() + self.compute_nplus()));
        assert!(nb_cells > 0);
        // store the nminus: nminus cells stay in the nminus
        // compartement, that is we dont apply Gaussian or Poisson
        // pdf to them
        let nminus = self.nminus;
        if self.compute_nplus() > 0 {
            match strategy {
                SamplingStrategy::Uniform => {}
                SamplingStrategy::Gaussian(_)
                | SamplingStrategy::Poisson(_)
                | SamplingStrategy::Exponential(_) => {
                    let mut ecdna_copies = Vec::with_capacity(
                        *self.get_nminus() as usize + self.compute_nplus() as usize,
                    );
                    // drop the nminus so the hist will not contain any nminus
                    // cells, the nminus cells stay in the 0 compartement (see
                    // comment above)
                    self.nminus = 0;
                    match strategy {
                        SamplingStrategy::Gaussian(sigma) => {
                            for (k, cells) in self.create_histogram() {
                                ecdna_copies.append(
                                    &mut rand_distr::Normal::new(k as f32, sigma.get())
                                        .unwrap()
                                        .sample_iter(rng.clone())
                                        .take(cells as usize)
                                        .map(|copy| f32::round(copy) as u16)
                                        .collect(),
                                );
                            }
                        }
                        SamplingStrategy::Poisson(rate) => {
                            for (k, cells) in self.create_histogram() {
                                ecdna_copies.append(
                                    &mut rand_distr::Poisson::new(k as f32 * rate.get())
                                        .unwrap()
                                        .sample_iter(rng.clone())
                                        .take(cells as usize)
                                        .map(|copy| f32::round(copy) as u16)
                                        .collect(),
                                );
                            }
                        }
                        SamplingStrategy::Exponential(rate) => {
                            let exp = rand_distr::Exp::new(rate.get()).unwrap();
                            for (k, cells) in self.create_histogram() {
                                ecdna_copies.append(
                                    &mut exp
                                        .sample_iter(rng.clone())
                                        .take(cells as usize)
                                        .map(|copy| f32::round(k as f32 * copy) as u16)
                                        .collect(),
                                );
                            }
                        }
                        _ => unreachable!(),
                    }

                    let distribution = EcDNADistribution::from(ecdna_copies);
                    // re-insert the nminus cells
                    self.nminus = nminus + distribution.nminus;
                    self.nplus = distribution.nplus;
                }
            }
            self.undersample(nb_cells, rng);
        } else {
            self.nminus = nb_cells;
        }
    }

    fn undersample(&mut self, nb_cells: u64, mut rng: impl Rng) {
        //! Draw a random sample without replacement from the
        //! `EcDNADistribution` by storing all cells into a `vec`, shuffling it
        //! and taking `nb_cells`.
        //!
        //! See [this](https://github.com/rust-random/book/blob/59649c93ed72e92c1644e3972a110f6ba5bc058d/src/guide-process.md).
        let mut distribution: Vec<u16> = self
            .nplus
            .iter()
            .map(|&k| k.get())
            .chain(std::iter::repeat(0u16).take(*self.get_nminus() as usize))
            .collect();
        let sample = if nb_cells as usize > self.nplus.len() / 2 {
            let tot = distribution.len();
            distribution
                .partial_shuffle(&mut rng, tot - nb_cells as usize)
                .1
        } else {
            distribution.partial_shuffle(&mut rng, nb_cells as usize).0
        };
        let new_distribution = EcDNADistribution::from(sample.to_vec());
        // re-shuffle because EcDNADistribution::from sort the nplus cells
        self.nplus = new_distribution.nplus;
        let amount = self.nplus.len() / 3;
        self.nplus.partial_shuffle(&mut rng, amount / 3);
        self.nminus = new_distribution.nminus;
    }

    pub fn ks_distance(&self, ecdna: &EcDNADistribution) -> f32 {
        //! The ks distance represents the maximal absolute distance between the
        //! empirical cumulative distributions of two `EcDNADistribution`s.
        //! ## Panics
        //! When the distributions are smaller than 8 samples.
        EcDNADistribution::calculate_statistic(&self.nplus, &ecdna.nplus)
    }

    fn calculate_statistic<T: Ord + Clone>(xs: &[T], ys: &[T]) -> f32 {
        // https://github.com/daithiocrualaoich/kolmogorov_smirnov/blob/cb067e92ec837efbad66e8bbcf85500ad778feb8/src/test.rs#L127
        assert!(!xs.is_empty());
        assert!(!ys.is_empty());
        let n = xs.len();
        let m = ys.len();

        assert!(n > 7 && m > 7);
        let mut xs = xs.to_vec();
        let mut ys = ys.to_vec();

        // xs and ys must be sorted for the stepwise ECDF calculations to work.
        xs.sort();
        ys.sort();

        // The current value testing for ECDF difference. Sweeps up through
        // elements present in xs and ys.
        let mut current: &T;

        // i, j index the first values in xs and ys that are greater than current.
        let mut i = 0;
        let mut j = 0;

        // ecdf_xs, ecdf_ys always hold the ECDF(current) of xs and ys.
        let mut ecdf_xs = 0.0;
        let mut ecdf_ys = 0.0;

        // The test statistic value computed over values <= current.
        let mut statistic = 0.0;

        while i < n && j < m {
            // Advance i through duplicate samples in xs.
            let x_i = &xs[i];

            while i + 1 < n && *x_i == xs[i + 1] {
                i += 1;
            }

            // Advance j through duplicate samples in ys.
            let y_j = &ys[j];

            while j + 1 < m && *y_j == ys[j + 1] {
                j += 1;
            }

            // Step to the next sample value in the ECDF sweep from low to high.
            current = min(x_i, y_j);

            // Update invariant conditions for i, j, ecdf_xs, and ecdf_ys.
            if current == x_i {
                ecdf_xs = (i + 1) as f32 / n as f32;
                i += 1;
            }

            if current == y_j {
                ecdf_ys = (j + 1) as f32 / m as f32;
                j += 1;
            }

            // Update invariant conditions for the test statistic.
            let diff = (ecdf_xs - ecdf_ys).abs();

            if diff > statistic {
                statistic = diff;
            }
        }

        // Don't need to walk the rest of the samples because one of the ecdfs is
        // already one and the other will be increasing up to one. This means the
        // difference will be monotonically decreasing, so we have our test
        // statistic value already.
        statistic
    }

    pub fn create_histogram(&self) -> HashMap<u16, u64> {
        let mut lookup = self.nplus.iter().fold(HashMap::new(), |mut acc, c| {
            let c_u16 = c.get();
            *acc.entry(c_u16).or_insert(0) += 1;
            acc
        });

        if self.nminus > 0 {
            lookup.insert(0, *self.get_nminus());
        }

        lookup
    }

    pub fn drop_cells_with_k_copies(self, k: DNACopy) -> EcDNADistribution {
        //! Create a new `EcDNADistribution` without cells with `k` copies.
        //!
        //! If you want to drop cells with 0 copies, use
        //! [`EcDNADistribution::drop_nminus`] which doesn't consume self.
        let nminus = self.nminus;
        let nplus = self.nplus.into_iter().filter(|&ecdna| ecdna != k).collect();
        EcDNADistribution { nplus, nminus }
    }

    pub fn drop_nminus(&mut self) {
        //! Drop all cells without any ecDNA copies (cells with 0 copies).
        self.nminus = 0;
    }

    pub fn scale_by(self, c: f32) -> Self {
        //! Scale (divide) the ecDNA distribution by a constant `c`.
        //!
        //! Since the ecDNA copies [`DNACopy`] are integers and, by dividing by
        //! a float, floats can arise, we round the scaled copies using
        //! ceil from [`std::f32`].
        //! ## Panics
        //! When `c` is smaller or equal than 0, is [`f32::NAN`] or
        //! [`f32::INFINITY`].
        if c.is_sign_positive() {
            assert!(!c.is_nan(), "`c` cannot be NaN");
            if c.is_finite() {
                if (c - 0f32).abs() < f32::EPSILON {
                    panic!("`c` cannot be zero!");
                } else {
                    // unchecked because we are using ceil and c cannot be inf.
                    let nplus = self
                        .nplus
                        .into_iter()
                        .map(|copy| unsafe {
                            NonZeroU16::new_unchecked((copy.get() as f32 / c).ceil() as u16)
                        })
                        .collect();
                    return EcDNADistribution {
                        nminus: self.nminus,
                        nplus,
                    };
                }
            }
            panic!("`c` cannot be infinite");
        }
        panic!("`c` must be greather than 0!")
    }

    pub fn compute_mean(&self) -> f32 {
        //! Returns `f32::NAN` if no cells are present.
        let mean = self.nplus.iter().map(|&ele| ele.get() as u32).sum::<u32>() as f32;
        let n = self.nminus as f32 + self.nplus.len() as f32;
        mean / n
    }

    pub fn compute_variance(&self) -> f32 {
        //! Compute the variance and returs `f32::NAN` if no cells are present.
        let mean = self.compute_mean();
        if mean.is_nan() {
            f32::NAN
        } else {
            let mut variance = self
                .nplus
                .iter()
                .map(|value| {
                    let diff = mean - (value.get() as f32);

                    diff * diff
                })
                .sum::<f32>();
            variance += *self.get_nminus() as f32 * mean * mean;
            variance /= (*self.get_nminus() + self.compute_nplus()) as f32;

            variance
        }
    }

    pub fn compute_entropy(&self) -> f32 {
        let ntot = self.compute_nplus() as f32 + self.nminus as f32;
        let lookup: HashMap<u16, f32> = self
            .create_histogram()
            .into_iter()
            .map(|(k, copy)| (k, copy as f32))
            .collect();

        -lookup
            .values()
            .map(|&count| {
                let prob = count / ntot;
                prob * prob.log2()
            })
            .sum::<f32>()
    }

    pub fn compute_frequency(&self) -> f32 {
        let nplus = self.compute_nplus() as f32;
        let nminus = self.nminus as f32;
        let freq = nplus / (nplus + nminus);
        if freq.is_nan() {
            0f32
        } else {
            freq
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::test_util::NonEmptyDistribtionWithNPlusCells;

    use super::*;
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;
    use std::{
        collections::{hash_map::RandomState, HashSet},
        num::{NonZeroU64, NonZeroU8},
    };

    #[derive(Clone, Debug)]
    pub struct NPlusVec(pub Vec<DNACopy>);

    impl Arbitrary for NPlusVec {
        fn arbitrary(g: &mut Gen) -> Self {
            NPlusVec(
                (0..4u16 * (u8::MAX as u16))
                    .map(|_| NonZeroU16::arbitrary(g))
                    .collect(),
            )
        }
    }

    #[derive(Clone, Debug)]
    pub struct SigmaGrZero(pub Sigma);

    impl Arbitrary for SigmaGrZero {
        fn arbitrary(g: &mut Gen) -> Self {
            SigmaGrZero(Sigma(LambdaGrZero::arbitrary(g).0))
        }
    }
    #[derive(Clone, Debug)]
    pub struct LambdaGrZero(pub Lambda);

    impl Arbitrary for LambdaGrZero {
        fn arbitrary(g: &mut Gen) -> Self {
            let mut lambda = u8::arbitrary(g) as f32 / 10.;
            while !lambda.is_normal() || lambda.is_sign_negative() {
                lambda = u8::arbitrary(g) as f32 / 10.;
            }
            LambdaGrZero(Lambda::new(lambda))
        }
    }

    #[test]
    fn ecdna_compute_mean_empty_test() {
        assert!(EcDNADistribution {
            nminus: 0,
            nplus: vec![]
        }
        .compute_mean()
        .is_nan());
    }

    #[test]
    fn ecdna_compute_mean_10() {
        assert_eq!(
            EcDNADistribution {
                nminus: 0,
                nplus: vec![
                    NonZeroU16::new(15).unwrap(),
                    NonZeroU16::new(5).unwrap(),
                    NonZeroU16::new(10).unwrap()
                ],
            }
            .compute_mean(),
            10f32
        )
    }

    #[test]
    fn ecdna_compute_mean_with_nminus_10() {
        assert_eq!(
            EcDNADistribution {
                nminus: 1,
                nplus: vec![
                    NonZeroU16::new(15).unwrap(),
                    NonZeroU16::new(15).unwrap(),
                    NonZeroU16::new(10).unwrap()
                ],
            }
            .compute_mean(),
            10f32
        )
    }

    #[test]
    fn ecdna_compute_mean_with_nminus_5() {
        assert_eq!(
            EcDNADistribution {
                nminus: 5,
                nplus: vec![
                    NonZeroU16::new(15).unwrap(),
                    NonZeroU16::new(15).unwrap(),
                    NonZeroU16::new(10).unwrap()
                ],
            }
            .compute_mean(),
            5f32
        )
    }

    #[test]
    fn ecdna_compute_mean_with_nminus_05() {
        assert_eq!(
            EcDNADistribution {
                nminus: 19,
                nplus: vec![NonZeroU16::new(10).unwrap()],
            }
            .compute_mean(),
            0.5f32
        )
    }

    #[quickcheck]
    fn ecdna_compute_mean_not_empty_test(copies: DNACopy, nminus: NonZeroU64) -> bool {
        let nminus = nminus.get();
        EcDNADistribution {
            nminus,
            nplus: vec![copies],
        }
        .compute_mean()
            == (copies.get() as f32 / nminus as f32)
    }

    #[quickcheck]
    fn ecdna_compute_mean_not_empty_no_nminus_test(copies: DNACopy) -> bool {
        EcDNADistribution {
            nminus: 0,
            nplus: vec![copies],
        }
        .compute_mean()
            == copies.get() as f32
    }

    #[test]
    fn compute_variance_no_cells() {
        let distribution = EcDNADistribution {
            nminus: 0,
            nplus: vec![],
        };
        assert!(distribution.compute_variance().is_nan());
    }

    #[quickcheck]
    fn compute_variance_no_nminus(nplus: NonZeroU8) -> bool {
        let distribution = EcDNADistribution {
            nminus: 0,
            nplus: vec![DNACopy::new(1).unwrap(); nplus.get() as usize],
        };
        (distribution.compute_variance() - 0.).abs() < f32::EPSILON
    }

    #[quickcheck]
    fn compute_variance_no_nplus(nminus: NonZeroU8) -> bool {
        let distribution = EcDNADistribution {
            nminus: nminus.get() as u64,
            nplus: vec![],
        };
        (distribution.compute_variance() - 0.).abs() < f32::EPSILON
    }

    #[test]
    fn compute_variance_equal_1() {
        let distribution = EcDNADistribution {
            nminus: 1,
            nplus: vec![DNACopy::new(2).unwrap()],
        };
        assert!((distribution.compute_variance() - 1.).abs() < f32::EPSILON);
    }

    #[quickcheck]
    fn ecdna_compute_frequency_no_nplus(nminus: u64) -> bool {
        (EcDNADistribution {
            nminus,
            nplus: vec![],
        }
        .compute_frequency()
            - 0f32)
            .abs()
            < f32::EPSILON
    }

    #[quickcheck]
    fn ecdna_compute_frequency_nonminus(copies: DNACopy) -> bool {
        (EcDNADistribution {
            nminus: 0,
            nplus: vec![copies],
        }
        .compute_frequency()
            - 1f32)
            .abs()
            < f32::EPSILON
    }

    #[quickcheck]
    fn ecdna_compute_frequency_01(copies: DNACopy) -> bool {
        (EcDNADistribution {
            nminus: 9,
            nplus: vec![copies],
        }
        .compute_frequency()
            - 0.1f32)
            .abs()
            < f32::EPSILON
    }

    #[quickcheck]
    fn ecdna_compute_frequency_02(copies: DNACopy) -> bool {
        (EcDNADistribution {
            nminus: 4,
            nplus: vec![copies],
        }
        .compute_frequency()
            - 0.2f32)
            .abs()
            < f32::EPSILON
    }

    #[quickcheck]
    fn ecdna_compute_frequency_05(copies: DNACopy) -> bool {
        (EcDNADistribution {
            nminus: 1,
            nplus: vec![copies],
        }
        .compute_frequency()
            - 0.5f32)
            .abs()
            < f32::EPSILON
    }

    #[test]
    fn compute_nplus_0_cells_test() {
        let ecdna = EcDNADistribution {
            nminus: 0,
            nplus: vec![],
        };
        assert_eq!(ecdna.compute_nplus(), 0);
        assert!(ecdna.is_empty());
    }

    #[quickcheck]
    fn compute_nplus_0_cells_with_nminus_test(nminus: NonZeroU8) -> bool {
        let ecdna = EcDNADistribution {
            nminus: nminus.get() as u64,
            nplus: vec![],
        };
        ecdna.compute_nplus() == 0 && !ecdna.is_empty()
    }

    #[quickcheck]
    fn compute_nplus_0_nminus_test(copies: DNACopy) -> bool {
        let ecdna = EcDNADistribution {
            nminus: 0,
            nplus: vec![copies],
        };
        ecdna.compute_nplus() > 0 && !ecdna.is_empty()
    }

    #[quickcheck]
    fn get_uniform_random_nplus_test(seed: u64, distr: NonEmptyDistribtionWithNPlusCells) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let idx = distr.0.pick_random_nplus(&mut rng);
        distr.0.nplus[idx].get() > 0
    }

    #[quickcheck]
    fn get_uniform_random_nplus_reproducible(
        seed: u64,
        distr: NonEmptyDistribtionWithNPlusCells,
    ) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let idx_first = distr.0.pick_random_nplus(&mut rng);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let idx_second = distr.0.pick_random_nplus(&mut rng);
        idx_first == idx_second && distr.0.nplus[idx_first].get() > 0
    }

    #[test]
    #[should_panic]
    fn ecdna_pick_random_nplus_panics_nminus_test() {
        let mut rng = ChaCha8Rng::from_entropy();
        let distr = EcDNADistribution {
            nminus: 12,
            nplus: vec![],
        };
        distr.pick_random_nplus(&mut rng);
    }

    #[quickcheck]
    fn pick_random_nplus_test(distribution: NonEmptyDistribtionWithNPlusCells, seed: u64) {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let random = distribution.0.pick_random_nplus(&mut rng);
        assert!(distribution.0.nplus.get(random).is_some());
    }

    #[quickcheck]
    fn pick_random_nplus_reproducible_test(
        distribution: NonEmptyDistribtionWithNPlusCells,
        seed: u64,
    ) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let idx_first = distribution.0.pick_random_nplus(&mut rng);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let idx_second = distribution.0.pick_random_nplus(&mut rng);

        idx_first == idx_second && distribution.0.nplus[idx_first].get() > 0
    }

    #[quickcheck]
    fn pick_random_nplus_no_choice_test(seed: u64, copies: DNACopy) -> bool {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let distr = EcDNADistribution {
            nminus: 12,
            nplus: vec![copies],
        };
        distr.nplus[distr.pick_random_nplus(&mut rng)] == copies
    }

    #[test]
    fn test_entropy_rosetta() {
        // https://rosettacode.org/wiki/Entropy#Rust
        let mut distribution = EcDNADistribution {
            nplus: vec![1, 2, 2, 3, 3, 3, 4, 4, 4, 4]
                .into_iter()
                .map(|ele| NonZeroU16::new(ele).unwrap())
                .collect(),
            nminus: 0,
        };

        let entropy = distribution.compute_entropy();
        let expected = 1.8464392f32;
        assert!(
            (entropy - expected).abs() < 0.0001,
            "should be {}, got {} instead",
            expected,
            entropy
        );

        distribution.nminus += 1;
        assert!(
            (distribution.compute_entropy() - 2.1180782).abs() < f32::EPSILON,
            "{}",
            entropy
        );
    }

    #[test]
    #[should_panic]
    fn ecdna_ks_distance_x_empty() {
        let x = EcDNADistribution {
            nplus: vec![],
            nminus: 0,
        };
        let y = EcDNADistribution {
            nplus: vec![NonZeroU16::new(3).unwrap(); 7],
            nminus: 0,
        };
        x.ks_distance(&y);
    }

    #[test]
    #[should_panic]
    fn ecdna_ks_distance_y_empty() {
        let x = EcDNADistribution {
            nplus: vec![NonZeroU16::new(3).unwrap(); 7],
            nminus: 0,
        };
        let y = EcDNADistribution {
            nplus: vec![],
            nminus: 0,
        };
        x.ks_distance(&y);
    }

    #[test]
    #[should_panic]
    fn ecdna_ks_distance_small_sample_y() {
        let x = EcDNADistribution {
            nplus: vec![
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(2).unwrap(),
            ],
            nminus: 1,
        };
        let y = EcDNADistribution {
            nplus: vec![NonZeroU16::new(1).unwrap(), NonZeroU16::new(2).unwrap()],
            nminus: 100,
        };
        x.ks_distance(&y);
    }

    #[test]
    #[should_panic]
    fn ecdna_ks_distance_small_sample_x() {
        let x = EcDNADistribution {
            nplus: vec![NonZeroU16::new(1).unwrap(), NonZeroU16::new(2).unwrap()],
            nminus: 100,
        };
        let y = EcDNADistribution {
            nplus: vec![
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(1).unwrap(),
                NonZeroU16::new(2).unwrap(),
            ],
            nminus: 100,
        };
        x.ks_distance(&y);
    }

    #[quickcheck]
    fn ecdna_ks_distance_same_data(x: NonEmptyDistribtionWithNPlusCells) -> bool {
        (x.0.ks_distance(&x.0) - 0f32).abs() <= f32::EPSILON
    }

    #[quickcheck]
    fn ecdna_ks_distance_is_one_when_no_overlap_in_support(
        x: NonEmptyDistribtionWithNPlusCells,
    ) -> bool {
        // https://github.com/daithiocrualaoich/kolmogorov_smirnov/blob/master/src/test.rs#L474
        let y_min = x.0.nplus.iter().max().unwrap().get() + 1;
        let y: EcDNADistribution = EcDNADistribution {
            nminus: 0,
            nplus: x
                .0
                .nplus
                .iter()
                .map(|ele| NonZeroU16::new(ele.get() + y_min).unwrap())
                .collect(),
        };
        (x.0.ks_distance(&y) - 1f32).abs() <= f32::EPSILON
    }

    #[quickcheck]
    fn ecdna_ks_distance_is_point_five_when_semi_overlap_in_support(
        x: NonEmptyDistribtionWithNPlusCells,
    ) -> bool {
        // https://github.com/daithiocrualaoich/kolmogorov_smirnov/blob/master/src/test.rs#L474
        let y_min = x.0.nplus.iter().max().unwrap().get() + 1;
        // Add all the original items back too.
        let mut y = EcDNADistribution {
            nminus: x.0.nminus,
            nplus: x.0.nplus.clone(),
        };

        for copy in x.0.nplus.iter() {
            y.nplus.push(NonZeroU16::new(copy.get() + y_min).unwrap());
        }
        y.nminus *= 2;
        let y_ntot = y.compute_nplus() + y.nminus;
        let ntot = x.0.compute_nplus() + x.0.nminus;
        assert_eq!(y_ntot, ntot * 2, "{} vs {}", y_ntot, ntot * 2);

        (x.0.ks_distance(&y) - 0.5f32).abs() <= f32::EPSILON
    }

    #[test]
    fn create_histogram_without_nminus() {
        let distribution = EcDNADistribution {
            nminus: 0,
            nplus: [1u16, 1u16, 2u16, 200u16]
                .iter()
                .map(|&ele| NonZeroU16::new(ele).unwrap())
                .collect(),
        };
        let map1: HashMap<u16, f32> = distribution
            .create_histogram()
            .into_iter()
            .map(|(k, copy)| (k, copy as f32))
            .collect();
        let map2 = HashMap::from([(1, 2f32), (2, 1f32), (200, 1f32)]);
        assert_eq!(map1.len(), map2.len());
        assert!(map1.keys().all(|k| map2.contains_key(k)));
        assert!(map1.keys().all(|k| map1.get(k) == map2.get(k)));
    }

    #[test]
    fn create_histogram() {
        let distribution = EcDNADistribution {
            nminus: 10,
            nplus: [1u16, 1u16, 2u16, 200u16]
                .iter()
                .map(|&ele| NonZeroU16::new(ele).unwrap())
                .collect(),
        };
        let map1: HashMap<u16, f32> = distribution
            .create_histogram()
            .into_iter()
            .map(|(k, copy)| (k, copy as f32))
            .collect();
        let map2 = HashMap::from([(0, 10f32), (1, 2f32), (2, 1f32), (200, 1f32)]);
        assert_eq!(map1.len(), map2.len());
        assert!(map1.keys().all(|k| map2.contains_key(k)));
        assert!(map1.keys().all(|k| map1.get(k) == map2.get(k)));
    }

    #[quickcheck]
    fn from_vec_empty_to_distribution(capacity: u8) -> bool {
        let distribution: EcDNADistribution = Vec::with_capacity(capacity as usize).into();
        distribution.nplus.capacity() == 0usize
            && distribution.nminus == 0u64
            && distribution.nplus.is_empty()
            && distribution.is_empty()
    }

    #[quickcheck]
    fn from_vec_to_distribution(mut nplus: NPlusVec, nminus: u8) -> bool {
        let mut my_vec: Vec<u16> = nplus.0.iter().map(|&k| k.get()).collect();
        my_vec.extend(std::iter::repeat(0u16).take(nminus as usize));
        my_vec.shrink_to_fit();
        let capacity = my_vec.capacity();
        let distribution: EcDNADistribution = my_vec.into();
        nplus.0.sort();
        distribution.nplus.capacity() <= capacity
            && distribution.nminus == nminus as u64
            && distribution.nplus == nplus.0
    }

    #[test]
    #[should_panic]
    fn sample_wrong_sample_size() {
        let my_vec = vec![0u16, 1u16];
        let mut distribution: EcDNADistribution = my_vec.into();
        distribution.sample(4, &SamplingStrategy::Uniform, ChaCha8Rng::seed_from_u64(26));
    }

    #[test]
    #[should_panic]
    fn sample_wrong_sample_size_0() {
        let my_vec = vec![0u16, 1u16];
        let mut distribution: EcDNADistribution = my_vec.into();
        distribution.sample(0, &SamplingStrategy::Uniform, ChaCha8Rng::seed_from_u64(26));
    }

    #[quickcheck]
    fn sample_gaussian_2_zeros_1_one_test(sigma: SigmaGrZero, seed: u64) -> bool {
        let distr_vec = vec![0, 0, 1];
        let size = distr_vec.len();
        let mut distribution = EcDNADistribution::from(distr_vec);
        distribution.sample(
            size as u64,
            &SamplingStrategy::Gaussian(sigma.0),
            ChaCha8Rng::seed_from_u64(seed),
        );
        distribution.compute_nplus() < 2
    }

    #[quickcheck]
    fn sample_gaussian_only_zeros_test(sigma: SigmaGrZero, seed: u64) -> bool {
        let distr_vec = vec![0; 100];
        let size = distr_vec.len();
        let mut distribution = EcDNADistribution::from(distr_vec);
        distribution.sample(
            size as u64,
            &SamplingStrategy::Gaussian(sigma.0),
            ChaCha8Rng::seed_from_u64(seed),
        );
        distribution.nplus.is_empty() && *distribution.get_nminus() as usize == size
    }

    #[quickcheck]
    fn sample_gaussian_range_copy_test(copy: DNACopy, seed: u64) -> bool {
        let range = 4;
        let mut min = 1;
        let max = copy.get() + range;
        if copy.get() > range {
            min = copy.get() - range;
        }
        let distr_vec = vec![copy.get(); 100];
        let size = distr_vec.len();
        let mut distribution = EcDNADistribution::from(distr_vec);
        distribution.sample(
            size as u64,
            &SamplingStrategy::Gaussian(Sigma::new(1.0)),
            ChaCha8Rng::seed_from_u64(seed),
        );
        distribution
            .nplus
            .into_iter()
            .all(|copy| min <= copy.get() && copy.get() <= max)
    }

    #[quickcheck]
    fn sample_gaussian_test(seed: u64, sigma: SigmaGrZero) -> bool {
        let distr_vec = vec![1; 100];
        let size = distr_vec.len();
        let mut distribution = EcDNADistribution::from(distr_vec);
        distribution.sample(
            size as u64,
            &SamplingStrategy::Gaussian(sigma.0),
            ChaCha8Rng::seed_from_u64(seed),
        );
        let mut found = false;
        for copy in distribution.nplus.iter() {
            if copy.get() >= 1 {
                found = true;
            }
        }
        found
    }

    #[should_panic]
    #[test]
    fn lambda_0_small_poisson_test() {
        Lambda::new(0.);
    }

    #[should_panic]
    #[test]
    fn lambda_too_small_poisson_test() {
        Lambda::new(-1.0);
    }

    #[quickcheck]
    fn sample_poisson_test(lambda: LambdaGrZero, seed: u64) -> bool {
        let distr_vec = vec![1; 100];
        let size = distr_vec.len();
        let mut distribution = EcDNADistribution::from(distr_vec);
        distribution.sample(
            size as u64,
            &SamplingStrategy::Poisson(lambda.0),
            ChaCha8Rng::seed_from_u64(seed),
        );
        distribution.get_nminus() + distribution.compute_nplus() == size as u64
    }

    #[quickcheck]
    fn sample_poisson_reproducible_test(
        lambda: LambdaGrZero,
        distribution: NonEmptyDistribtionWithNPlusCells,
        seed: u64,
    ) -> bool {
        let size = distribution.0.compute_nplus() + *distribution.0.get_nminus();

        let mut distr1 = distribution.0.clone();
        distr1.sample(
            size,
            &SamplingStrategy::Poisson(lambda.0),
            ChaCha8Rng::seed_from_u64(seed),
        );

        let mut distr2 = distribution.0;
        distr2.sample(
            size,
            &SamplingStrategy::Poisson(lambda.0),
            ChaCha8Rng::seed_from_u64(seed),
        );

        distr1 == distr2
    }

    #[should_panic]
    #[test]
    fn lambda_0_small_exp_test() {
        Lambda::new(0.);
    }

    #[should_panic]
    #[test]
    fn lambda_too_small_exp_test() {
        Lambda::new(-1.0);
    }

    #[quickcheck]
    fn sample_exp_test(lambda: LambdaGrZero, seed: u64) -> bool {
        let distr_vec = vec![1; 100];
        let size = distr_vec.len();
        let mut distribution = EcDNADistribution::from(distr_vec);
        distribution.sample(
            size as u64,
            &SamplingStrategy::Exponential(lambda.0),
            ChaCha8Rng::seed_from_u64(seed),
        );
        distribution.get_nminus() + distribution.compute_nplus() == size as u64
    }

    #[derive(Debug, Clone)]
    enum SamplingStrategyEnum {
        Gaussian,
        Poisson,
        Exp,
    }

    impl Arbitrary for SamplingStrategyEnum {
        fn arbitrary(g: &mut Gen) -> Self {
            match u8::arbitrary(g) % 3 {
                0 => Self::Gaussian,
                1 => Self::Poisson,
                2 => Self::Exp,
                _ => panic!(),
            }
        }
    }

    #[quickcheck]
    fn sample_reproducible_test(
        parameter: LambdaGrZero,
        sampling_strategy: SamplingStrategyEnum,
        distribution: NonEmptyDistribtionWithNPlusCells,
        seed: u64,
    ) -> bool {
        let size = distribution.0.compute_nplus() + *distribution.0.get_nminus();
        let sampling = match sampling_strategy {
            SamplingStrategyEnum::Gaussian => {
                SamplingStrategy::Gaussian(Sigma::new(parameter.0.get()))
            }
            SamplingStrategyEnum::Poisson => SamplingStrategy::Poisson(parameter.0),
            SamplingStrategyEnum::Exp => SamplingStrategy::Exponential(parameter.0),
        };

        let mut distr1 = distribution.0.clone();
        distr1.sample(size, &sampling, ChaCha8Rng::seed_from_u64(seed));

        let mut distr2 = distribution.0;
        distr2.sample(size, &sampling, ChaCha8Rng::seed_from_u64(seed));

        distr1 == distr2
    }

    #[quickcheck]
    fn undersample_full_distribution(
        mut distribution: NonEmptyDistribtionWithNPlusCells,
        seed: u64,
    ) -> bool {
        let size = (distribution.0.compute_nplus() + *distribution.0.get_nminus()) as usize;
        let mut distribution_sampled = distribution.0.clone();
        distribution_sampled.undersample(size as u64, &mut ChaCha8Rng::seed_from_u64(seed));
        distribution.0.nplus.sort();
        distribution_sampled.nplus.sort();

        distribution.0 == distribution_sampled
    }

    #[quickcheck]
    fn undersample_distribution(
        size: NonZeroU8,
        mut distribution: NonEmptyDistribtionWithNPlusCells,
        seed: u64,
    ) -> bool {
        let copies_present = HashSet::<u16, RandomState>::from_iter(
            distribution
                .0
                .nplus
                .iter()
                .map(|&ele| ele.get())
                .chain(std::iter::repeat(0u16).take(size.get() as usize)),
        );

        distribution
            .0
            .undersample(size.get() as u64, &mut ChaCha8Rng::seed_from_u64(seed));
        let expected_size = size.get() as u64;
        let distr_size = *distribution.0.get_nminus() + distribution.0.compute_nplus();

        expected_size == distr_size
            && distribution
                .0
                .nplus
                .iter()
                .all(|&ele| copies_present.contains(&ele.get()))
    }

    #[test]
    #[should_panic]
    fn scale_by_zero_test() {
        let distribution = EcDNADistribution::from(vec![0, 1, 1, 2, 3, 4]);
        distribution.scale_by(0f32);
    }

    #[test]
    #[should_panic]
    fn scale_by_nan_test() {
        let distribution = EcDNADistribution::from(vec![0, 1, 1, 2, 3, 4]);
        distribution.scale_by(f32::NAN);
    }

    #[test]
    #[should_panic]
    fn scale_by_inf_test() {
        let distribution = EcDNADistribution::from(vec![0, 1, 1, 2, 3, 4]);
        distribution.scale_by(f32::INFINITY);
    }

    #[test]
    #[should_panic]
    fn scale_by_minus_inf_test() {
        let distribution = EcDNADistribution::from(vec![0, 1, 1, 2, 3, 4]);
        distribution.scale_by(f32::NEG_INFINITY);
    }

    #[test]
    #[should_panic]
    fn scale_by_neg_test() {
        let distribution = EcDNADistribution::from(vec![0, 1, 1, 2, 3, 4]);
        distribution.scale_by(-0f32);
    }

    #[test]
    fn scale_by_c_greather_than_1_test() {
        let mut distribution = EcDNADistribution::from(vec![0, 1, 1, 2, 3, 4]);
        distribution = distribution.scale_by(2f32);
        let expected = EcDNADistribution::from(vec![0, 1, 1, 1, 2, 2]);
        assert_eq!(distribution, expected);
    }

    #[test]
    fn scale_by_c_smaller_than_1_test() {
        let mut distribution = EcDNADistribution::from(vec![0, 1, 1, 2, 3, 4]);
        distribution = distribution.scale_by(0.5f32);
        let expected = EcDNADistribution::from(vec![0, 2, 2, 4, 6, 8]);
        assert_eq!(distribution, expected);
    }

    #[quickcheck]
    fn test_drop_cells_with_k_copies(distribution: NonEmptyDistribtionWithNPlusCells) -> bool {
        let k = distribution
            .0
            .nplus
            .choose(&mut rand::thread_rng())
            .unwrap()
            .to_owned();
        let filtered_distr = distribution.0.drop_cells_with_k_copies(k);

        for cell in filtered_distr.nplus {
            if cell == k {
                return false;
            }
        }
        true
    }

    #[quickcheck]
    fn test_display(distribution: NonEmptyDistribtionWithNPlusCells) {
        println!("{}", distribution.0);
    }
}
