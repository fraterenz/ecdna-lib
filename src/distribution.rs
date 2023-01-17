use std::{collections::HashMap, fs, num::NonZeroU16, path::Path};

use anyhow::{bail, Context};
use rand::{seq::SliceRandom, Rng};
use rand_chacha::ChaCha8Rng;

use crate::DNACopy;

const EPSILON_KS_STAT: f32 = f32::EPSILON * 10f32;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct EcDNADistribution {
    nminus: u64,
    nplus: Vec<DNACopy>,
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

    pub fn pick_remove_random_nplus(&mut self, rng: &mut ChaCha8Rng) -> anyhow::Result<DNACopy> {
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

    fn pick_random_nplus(&self, rng: &mut ChaCha8Rng) -> usize {
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

    pub fn decrease_nplus(&mut self, rng: &mut ChaCha8Rng, verbosity: u8) {
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
        let map = self.create_histogram();
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

    pub fn undersample(&mut self, nb_cells: u64, rng: &mut ChaCha8Rng) {
        //! Draw a random sample without replacement from the
        //! `EcDNADistribution` by storing all cells into a `vec`, shuffling it
        //! and taking `nb_cells`.
        //!
        //! See [this](https://github.com/rust-random/book/blob/59649c93ed72e92c1644e3972a110f6ba5bc058d/src/guide-process.md).
        //!
        //! ## Panics
        //! Panics when `nb_cells` is bigger or equal to the cells in the
        //! distribution or wheb `nb_cells` is 0.
        assert!(nb_cells < (*self.get_nminus() + self.compute_nplus()));
        assert!(nb_cells > 0);
        let mut distribution: Vec<u16> = self
            .nplus
            .iter()
            .map(|&k| k.get())
            .chain(std::iter::repeat(0u16).take(*self.get_nminus() as usize))
            .collect();
        // shuffle and take the first `nb_cells`
        let (sample, _) = distribution.partial_shuffle(rng, nb_cells as usize);
        let new_distribution = EcDNADistribution::from(sample.to_vec());
        self.nplus = new_distribution.nplus;
        self.nminus = new_distribution.nminus;
    }

    pub fn ks_distance(&self, ecdna: &EcDNADistribution) -> (f32, bool) {
        //! The ks distance represents the maximal absolute distance between the
        //! empirical cumulative distributions of two `EcDNADistribution`s.
        //!
        //! Compute ks distance with `NPlus` and `NMinus` cells, and returns the
        //! distance as well whether the loop stopped before reaching the maximal
        //! allowed copy number `u16::MAX`, which is the case when the distance is
        //! 1 ie maximal, or one of the cumulative distributions have reached 1
        //! and thus the distance can only decrease monotonically.
        //!
        //! Does **not** panic if empty distributions.
        // do not test small samples because ks is not reliable (ecdf)
        let ecdna_ntot = ecdna.compute_nplus() + ecdna.nminus;
        let ntot = self.compute_nplus() + self.nminus;
        if ntot < 10 || ecdna_ntot < 10 {
            return (1f32, false);
        }

        let mut distance = 0f32;
        // Compare the two empirical cumulative distributions (self) and ecdna
        let mut ecdf = 0f32;
        let mut ecdf_other = 0f32;
        let distribution = self.create_histogram_f32();
        let distribution_other = ecdna.create_histogram_f32();

        // iter over all ecDNA copies present in both data (self) and ecdna
        for copy in 0u16..u16::MAX {
            if let Some(ecdna_copy) = distribution.get(&copy) {
                ecdf += ecdna_copy / (ntot as f32);
            }
            if let Some(ecdna_copy) = distribution_other.get(&copy) {
                ecdf_other += ecdna_copy / (ecdna_ntot as f32);
            }
            // store the maximal distance between the two ecdfs
            let diff = (ecdf - ecdf_other).abs();
            if diff - distance >= EPSILON_KS_STAT {
                distance = diff;
            }
            // check if any of the ecdf have reached 1. If it's the case
            // the difference will decrease monotonically and we can stop
            let max_dist = (ecdf - 1.0).abs() <= EPSILON_KS_STAT
                || (ecdf_other - 1.0).abs() <= EPSILON_KS_STAT
                || (distance - 1.0).abs() <= EPSILON_KS_STAT;

            if max_dist {
                return (distance, true);
            }
        }

        (distance, false)
    }

    fn create_histogram(&self) -> HashMap<u16, u64> {
        let mut lookup = self.nplus.iter().fold(HashMap::new(), |mut acc, c| {
            let c_u16 = c.get();
            *acc.entry(c_u16).or_insert(0) += 1u64;
            acc
        });

        if self.nminus > 0 {
            lookup.insert(0, self.nminus);
        }

        lookup
    }

    fn create_histogram_f32(&self) -> HashMap<u16, f32> {
        let mut lookup = self.nplus.iter().fold(HashMap::new(), |mut acc, c| {
            let c_u16 = c.get();
            *acc.entry(c_u16).or_insert(0f32) += 1f32;
            acc
        });

        if self.nminus > 0 {
            lookup.insert(0, self.nminus as f32);
        }
        lookup
    }

    pub fn drop_nminus(&mut self) {
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

    pub fn compute_entropy(&self) -> f32 {
        let ntot = self.compute_nplus() as f32 + self.nminus as f32;
        let lookup = self.create_histogram_f32();

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
        assert!(ecdna.compute_nplus() == 0 && ecdna.is_empty());
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
    fn ecdna_ks_distance_empty() {
        let x = EcDNADistribution {
            nplus: vec![],
            nminus: 0,
        };
        let y = EcDNADistribution {
            nplus: vec![],
            nminus: 0,
        };
        assert_eq!(x.ks_distance(&y), (1f32, false));
        assert_eq!(y.ks_distance(&x), (1f32, false));
        assert_eq!(x.ks_distance(&x), (1f32, false));
    }

    #[test]

    fn ecdna_ks_distance_small_samples() {
        let x = EcDNADistribution {
            nplus: vec![NonZeroU16::new(1).unwrap(), NonZeroU16::new(2).unwrap()],
            nminus: 4,
        };
        let y = EcDNADistribution {
            nplus: vec![NonZeroU16::new(1).unwrap(), NonZeroU16::new(2).unwrap()],
            nminus: 3,
        };
        assert_eq!(x.ks_distance(&y), (1f32, false));
        assert_eq!(y.ks_distance(&x), (1f32, false));
        assert_eq!(x.ks_distance(&x), (1f32, false));
    }

    #[quickcheck]
    fn ecdna_ks_distance_same_data(x: NonEmptyDistribtionWithNPlusCells) -> bool {
        let (distance, _) = x.0.ks_distance(&x.0);
        (distance - 0f32).abs() <= f32::EPSILON
    }

    #[quickcheck]
    // #[ignore = "TODO"]
    fn ecdna_ks_distance_max_copy_number_u16(
        mut x: NonEmptyDistribtionWithNPlusCells,
        y: NonEmptyDistribtionWithNPlusCells,
    ) -> bool {
        x.0.nplus.push(NonZeroU16::new(u16::MAX).unwrap());
        let (_, convergence) = x.0.ks_distance(&y.0);
        convergence
    }

    #[quickcheck]
    // #[ignore = "Small error: not exacty 1"]
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
        let (distance, convergence) = x.0.ks_distance(&y);
        assert!(x.0.nminus > 0);
        (distance - 1f32).abs() <= EPSILON_KS_STAT && convergence
    }

    #[quickcheck]
    #[ignore = "Small error: not exacty 0.5"]
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

        let (distance, convergence) = x.0.ks_distance(&y);
        (distance - 0.5f32).abs() <= EPSILON_KS_STAT && convergence
    }

    #[quickcheck]
    fn ecdna_ks_distance_fast_convergence(
        x: NonEmptyDistribtionWithNPlusCells,
        nminus: NonZeroU8,
    ) -> bool {
        let y = EcDNADistribution {
            nminus: nminus.get() as u64,
            nplus: vec![NonZeroU16::new(1).unwrap(); 100],
        };
        let (_, convergence) = x.0.ks_distance(&y);
        let (_, convergence1) = y.ks_distance(&x.0);
        assert!(x.0.nplus.len() >= 10, "{} len", x.0.nplus.len());
        convergence && convergence1
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
        let map1 = distribution.create_histogram_f32();
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
        let map1 = distribution.create_histogram_f32();
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
    fn undersample_wrong_sample_size() {
        let my_vec = vec![0u16, 1u16];
        let mut distribution: EcDNADistribution = my_vec.into();
        distribution.undersample(4, &mut ChaCha8Rng::seed_from_u64(26));
    }

    #[test]
    #[should_panic]
    fn undersample_wrong_sample_size_0() {
        let my_vec = vec![0u16, 1u16];
        let mut distribution: EcDNADistribution = my_vec.into();
        distribution.undersample(0, &mut ChaCha8Rng::seed_from_u64(26));
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
}
