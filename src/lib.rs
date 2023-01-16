use std::num::NonZeroU16;

/// The ecDNA distribution is a collection of cells with and without ecDNA
/// copies.
pub mod distribution;

/// Approximate Bayesian inference performed on two [`distribution::EcDNADistribution`]s.
pub mod abc;

/// EcDNA copies are by definition non-zero.
/// We assume that the maximal ecDNA copies present in a system cannot be
/// greater than 65535 copies (u16 is 2^16 - 1 = 65535 copies).
pub type DNACopy = NonZeroU16;

#[cfg(test)]
pub mod test_util {
    use std::{collections::HashMap, num::NonZeroU8};

    use quickcheck::{Arbitrary, Gen};

    use crate::{distribution::EcDNADistribution, DNACopy};

    #[derive(Clone, Debug)]
    pub struct DNACopySegregatingGreatherThanOne(pub DNACopy);

    impl Arbitrary for DNACopySegregatingGreatherThanOne {
        fn arbitrary(g: &mut Gen) -> DNACopySegregatingGreatherThanOne {
            let mut copy = DNACopy::new(NonZeroU8::arbitrary(g).get() as u16).unwrap();
            if copy == DNACopy::new(1).unwrap() {
                copy = DNACopy::new(2).unwrap();
            }
            DNACopySegregatingGreatherThanOne(copy)
        }
    }

    #[derive(Clone, Debug)]
    pub struct NonEmptyDistribtionWithNPlusCells(pub EcDNADistribution);

    impl Arbitrary for NonEmptyDistribtionWithNPlusCells {
        fn arbitrary(g: &mut Gen) -> NonEmptyDistribtionWithNPlusCells {
            const MAX_ENTRIES: usize = 500;
            let mut distr = HashMap::with_capacity(MAX_ENTRIES);
            for _ in 0..MAX_ENTRIES {
                let copy = DNACopySegregatingGreatherThanOne::arbitrary(g);
                let cells = NonZeroU8::arbitrary(g).get() as u64;
                distr.insert(u16::from(copy.0), cells);
            }
            let cells = NonZeroU8::arbitrary(g).get() as u64;
            distr.insert(0, cells);
            let distr = EcDNADistribution::new(distr, 1000);
            NonEmptyDistribtionWithNPlusCells(distr)
        }
    }
}
