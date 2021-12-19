use ark_bls12_381::Fr;
use ark_ff::UniformRand;
use ark_serialize::CanonicalSerialize;
use blake2::{Blake2b, Digest};
use kgz::KzgCommitment;
use rand::{prelude::StdRng, SeedableRng};

#[derive(Clone)]
pub struct ChallengeGenerator {
    data: Vec<u8>,
}

impl ChallengeGenerator {
    pub fn new() -> Self {
        Self { data: vec![] }
    }
    pub fn digest(&mut self, commitment: &KzgCommitment) {
        commitment
            .inner()
            .serialize_unchecked(&mut self.data)
            .unwrap();
    }
    pub fn with_digest(commitments: &[KzgCommitment]) -> Self {
        let mut generator = Self::new();
        for commitment in commitments {
            generator.digest(commitment);
        }
        generator
    }
    fn generate_rng(self) -> StdRng {
        let mut hasher = Blake2b::new();
        let Self { data } = self;
        hasher.update(data);
        let hash: Vec<u8> = hasher.finalize().to_vec();
        let mut seed: [u8; 8] = Default::default();
        seed.copy_from_slice(&hash[0..8]);
        let seed = u64::from_le_bytes(seed);
        StdRng::seed_from_u64(seed)
    }
    pub fn generate_challenges<const N: usize>(self) -> [Fr; N] {
        let mut rng = self.generate_rng();

        let points = [0; N];
        points.map(|_| Fr::rand(&mut rng))
    }
}
