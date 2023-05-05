use anyhow::{ensure, Context};
use hdrhistogram::Histogram;
use serde::{ser::SerializeStruct, Deserialize, Serialize, Serializer};
use std::{fs, path::Path};

use derive_builder::Builder;

use crate::distribution::EcDNADistribution;

#[derive(Debug, Deserialize)]
pub struct ABCResultFitness {
    pub result: ABCResult,
    pub rates: Vec<f32>,
}

impl Serialize for ABCResultFitness {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("ABCResultFitness", 12)?;
        state.serialize_field("idx", &self.result.idx)?;
        state.serialize_field("mean", &self.result.mean)?;
        state.serialize_field("mean_stat", &self.result.mean_stat)?;
        state.serialize_field("frequency", &self.result.frequency)?;
        state.serialize_field("frequency_stat", &self.result.frequency_stat)?;
        state.serialize_field("entropy", &self.result.entropy)?;
        state.serialize_field("entropy_stat", &self.result.entropy_stat)?;
        state.serialize_field("ecdna_stat", &self.result.ecdna_stat)?;
        state.serialize_field("pop_size", &self.result.pop_size)?;
        state.serialize_field("b0", &self.rates[0])?;
        state.serialize_field("b1", &self.rates[1])?;
        if self.rates.len() > 2 {
            state.serialize_field("d0", &self.rates[2])?;
            state.serialize_field("d1", &self.rates[3])?;
            if self.rates.len() > 4 {
                unreachable!()
            }
        } else {
            state.serialize_field("d0", &0f32)?;
            state.serialize_field("d1", &0f32)?;
        }
        state.serialize_field("dropped_nminus", &self.result.dropped_nminus)?;

        state.end()
    }
}

pub struct ABCResultsFitness(pub Vec<ABCResultFitness>);

impl ABCResultsFitness {
    pub fn save(self, path2folder: &Path, verbosity: u8) -> anyhow::Result<()> {
        fs::create_dir_all(path2folder).with_context(|| "Cannot create dir".to_string())?;

        let mut abc = path2folder.join("abc");
        abc.set_extension("csv");
        if verbosity > 1 {
            println!("Saving ABC results to {:#?}", abc);
        }
        let mut wtr = csv::Writer::from_path(abc)?;

        for res in self.0 {
            if verbosity > 0 {
                println!(
                    "Statistics of run {}: Mean: {}, Freq: {}, Entropy: {}",
                    res.result.idx, res.result.mean, res.result.frequency, res.result.entropy
                );
            }
            wtr.serialize(&res)
                .with_context(|| "Cannot serialize the results from ABC inference".to_string())?;
        }
        wtr.flush()?;
        Ok(())
    }

    pub fn from_csv(path2csv: &Path) -> anyhow::Result<Self> {
        ensure!(path2csv.is_file());
        ensure!(path2csv.extension().unwrap() == "csv");

        let mut records: Vec<ABCResultFitness> = Vec::new();

        let mut rdr = csv::Reader::from_path(path2csv)
            .with_context(|| format!("cannot open csv file {:#?}", path2csv))?;
        for result in rdr.deserialize() {
            let record = result?;
            records.push(record);
        }
        Ok(ABCResultsFitness(records))
    }
}

/// The data used to perform ABC.
#[derive(Debug)]
pub struct Data {
    pub distribution: Option<EcDNADistribution>,
    pub mean: Option<f32>,
    pub frequency: Option<f32>,
    pub entropy: Option<f32>,
}

/// Perform the ABC rejection algorithm for one run to infer the most probable
/// values of the rates based on the patient's data.
///
/// ABC infer the most probable values of the birth-death rates (proliferation
/// rates and death rates) by comparing the summary statistics of the run
/// generated by the birth-death process against the summary statistics of the
/// patient's data.
///
/// When testing multiple statistics with ABC, save runs only if all statistics
/// pass the tests
pub struct ABCRejection;

impl ABCRejection {
    pub fn run(
        mut builder: ABCResultBuilder,
        distribution: &EcDNADistribution,
        target: &Data,
        drop_nminus: bool,
        verbosity: u8,
    ) -> ABCResult {
        //! Run the ABC rejection method by comparing a ecDNA `distribution`
        //! against the `target`.
        let pop_size = distribution.get_nminus() + distribution.compute_nplus();
        builder.pop_size(pop_size);
        if let Some(target_distribution) = &target.distribution {
            if target_distribution.compute_nplus() <= 7 || distribution.compute_nplus() <= 7 {
                if verbosity > 0 {
                    println!(
                        "Not enough samples in the distributions: {} and {}",
                        target_distribution.compute_nplus(),
                        distribution.compute_nplus()
                    );
                }
            } else {
                builder.ecdna_stat(Some(target_distribution.ks_distance(distribution)));
            }
        };

        let frequency = distribution.compute_frequency();
        builder.frequency(frequency);
        let frequency_stat = target.frequency.as_ref().map(|target_frequency| {
            if !drop_nminus {
                relative_change(target_frequency, &frequency)
            } else {
                f32::NAN
            }
        });
        builder.frequency_stat(frequency_stat);

        let entropy = distribution.compute_entropy();
        builder.entropy(entropy);
        let entropy_stat = target
            .entropy
            .as_ref()
            .map(|target_entropy| relative_change(target_entropy, &entropy));
        builder.entropy_stat(entropy_stat);

        let k_upper_bound = 4000;
        let mut hist = Histogram::<u64>::new_with_max(k_upper_bound, 2).unwrap();
        for (ecdna, nb_cells) in distribution.create_histogram().into_iter() {
            hist.record_n(ecdna as u64, nb_cells)
                .expect("should be in range");
        }

        let mean = hist.mean() as f32;
        builder.mean(mean);
        let mean_stat = target
            .mean
            .as_ref()
            .map(|target_mean| relative_change(target_mean, &mean));
        builder.mean_stat(mean_stat);

        if let Some(target_distribution) = target.distribution.as_ref() {
            let quantile = 0.9;

            let mut hist_target = Histogram::<u64>::new_with_max(k_upper_bound, 2).unwrap();
            for (ecdna, nb_cells) in target_distribution.create_histogram().into_iter() {
                hist_target
                    .record_n(ecdna as u64, nb_cells)
                    .expect("should be in range");
            }

            builder.k_max(Some(relative_change(
                &(hist_target.value_at_quantile(quantile) as f32),
                &(hist.value_at_quantile(quantile) as f32),
            )));
        }

        if verbosity > 0 {
            println!(
                "The stats are: ks:{:#?}, mean: {:#?}, freq: {:#?}, entropy: {:#?}, k_max: {:#?}",
                builder.ecdna_stat.unwrap_or(None),
                mean_stat,
                frequency_stat,
                entropy_stat,
                builder.k_max.unwrap_or(None)
            );
        }
        builder.dropped_nminus(drop_nminus);

        builder.build().expect("Cannot build ABC results")
    }
}

#[derive(Builder, Debug, Deserialize)]
pub struct ABCResult {
    pub idx: usize,
    pub mean: f32,
    #[builder(default)]
    pub mean_stat: Option<f32>,
    pub frequency: f32,
    #[builder(default)]
    pub frequency_stat: Option<f32>,
    pub entropy: f32,
    #[builder(default)]
    pub entropy_stat: Option<f32>,
    #[builder(default)]
    pub ecdna_stat: Option<f32>,
    #[builder(default)]
    pub k_max: Option<f32>,
    pub pop_size: u64,
    pub dropped_nminus: bool,
}

/// Relative change between two scalars
pub fn relative_change(x1: &f32, &x2: &f32) -> f32 {
    (x1 - x2).abs() / x1
}

#[cfg(test)]
mod tests {
    use quickcheck_macros::quickcheck;

    use crate::test_util::NonEmptyDistribtionWithNPlusCells;

    use super::*;

    #[quickcheck]
    fn abc_run_small_sample_with_nminus_target(
        distribution: NonEmptyDistribtionWithNPlusCells,
        drop_nminus: bool,
    ) -> bool {
        let target_distr = EcDNADistribution::from(vec![1, 2, 0]);
        let mut builder = ABCResultBuilder::default();
        builder.idx(1);
        let target = Data {
            distribution: Some(target_distr),
            mean: None,
            frequency: None,
            entropy: None,
        };
        let results = ABCRejection::run(builder, &distribution.0, &target, drop_nminus, 0);
        results.ecdna_stat.is_none() && results.dropped_nminus == drop_nminus
    }

    #[quickcheck]
    fn abc_run_small_sample_with_nminus(
        target_distr: NonEmptyDistribtionWithNPlusCells,
        drop_nminus: bool,
    ) -> bool {
        let distribution = EcDNADistribution::from(vec![1, 2, 0]);
        let mut builder = ABCResultBuilder::default();
        builder.idx(1);
        let target = Data {
            distribution: Some(target_distr.0),
            mean: None,
            frequency: None,
            entropy: None,
        };
        let results = ABCRejection::run(builder, &distribution, &target, drop_nminus, 0);
        results.ecdna_stat.is_none() && results.dropped_nminus == drop_nminus
    }

    #[quickcheck]
    fn abc_run_test(
        distribution: NonEmptyDistribtionWithNPlusCells,
        idx: usize,
        drop_nminus: bool,
    ) -> bool {
        let mut builder = ABCResultBuilder::default();
        builder.idx(idx);
        let mean = distribution.0.compute_mean();
        let frequency = distribution.0.compute_frequency();
        let entropy = distribution.0.compute_entropy();
        builder.mean(mean);
        builder.frequency(frequency);
        builder.entropy(entropy);
        let target = Data {
            distribution: Some(distribution.0.clone()),
            mean: Some(mean),
            frequency: Some(frequency),
            entropy: Some(entropy),
        };

        let results = ABCRejection::run(builder, &distribution.0, &target, drop_nminus, 0);
        let freq_test = if !drop_nminus {
            (results.frequency_stat.unwrap() - 0f32).abs() < f32::EPSILON
        } else {
            results.frequency_stat.unwrap().is_nan()
        };
        dbg!(&results.k_max);
        (results.ecdna_stat.unwrap() - 0f32).abs() < f32::EPSILON
            && (results.mean_stat.unwrap() - 0f32).abs() < f32::EPSILON
            && freq_test
            && (results.entropy_stat.unwrap() - 0f32).abs() < 10f32 * f32::EPSILON
            && results.dropped_nminus == drop_nminus
            && (results.k_max.unwrap() - 0f32).abs() < f32::EPSILON
    }

    #[quickcheck]
    fn abc_run_no_distribution_test(
        distribution: NonEmptyDistribtionWithNPlusCells,
        idx: usize,
        drop_nminus: bool,
    ) -> bool {
        let mut builder = ABCResultBuilder::default();
        builder.idx(idx);
        let mean = distribution.0.compute_mean();
        let frequency = distribution.0.compute_frequency();
        let entropy = distribution.0.compute_entropy();
        builder.mean(mean);
        builder.frequency(frequency);
        builder.entropy(entropy);
        let target = Data {
            distribution: None,
            mean: Some(mean),
            frequency: Some(frequency),
            entropy: Some(entropy),
        };

        let results = ABCRejection::run(builder, &distribution.0, &target, drop_nminus, 0);
        let freq_test = if !drop_nminus {
            (results.frequency_stat.unwrap() - 0f32).abs() < f32::EPSILON
        } else {
            results.frequency_stat.unwrap().is_nan()
        };
        results.ecdna_stat.is_none()
            && (results.mean_stat.unwrap() - 0f32).abs() < f32::EPSILON
            && freq_test
            && (results.entropy_stat.unwrap() - 0f32).abs() < 10f32 * f32::EPSILON
            && results.dropped_nminus == drop_nminus
            && results.k_max.is_none()
    }

    #[quickcheck]
    fn abc_run_distribution_only_test(
        distribution: NonEmptyDistribtionWithNPlusCells,
        idx: usize,
        drop_nminus: bool,
    ) -> bool {
        let mut builder = ABCResultBuilder::default();
        builder.idx(idx);
        let target = Data {
            distribution: Some(distribution.0.clone()),
            mean: None,
            frequency: None,
            entropy: None,
        };

        let results = ABCRejection::run(builder, &distribution.0, &target, drop_nminus, 0);
        (results.ecdna_stat.unwrap() - 0f32).abs() < f32::EPSILON
            && results.mean_stat.is_none()
            && results.frequency_stat.is_none()
            && results.entropy_stat.is_none()
            && results.dropped_nminus == drop_nminus
    }
}
