use crate::{Circuit, Poly, Prof};
use ark_bls12_381::Fr;
use ark_ff::{bytes::ToBytes, BigInteger256, One, UniformRand};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use blake2::{Blake2b, Digest};
use kgz::{KzgCommitment, KzgOpening, KzgScheme};
use rand::{prelude::StdRng, Rng, SeedableRng};
use std::{convert::TryInto, io::Write, iter::repeat_with};

use self::challenges::ChallengeGenerator;
mod challenges;

impl Circuit {
    fn proof(&self, a: Poly, b: Poly, c: Poly) -> Prof {
        todo!()
    }
}

struct PolyProof {
    commitment: KzgCommitment,
    opening: KzgOpening,
    evaluation: Fr,
}
struct PermutationProof {
    commitment: KzgCommitment,
    first_opening: KzgOpening,
    first_evaluation: Fr,
    second_opening: KzgOpening,
    second_evaluation: Fr,
}
struct Proof {
    a: PolyProof,
    b: PolyProof,
    c: PolyProof,
    permutation: PermutationProof,
}
fn prove(circuit: &Circuit, advice: [Poly; 3]) -> Proof {
    let scheme = KzgScheme::new(&circuit.srs);
    let domain = &circuit.domain;
    let rows = circuit.rows;

    let [a, b, c] = advice;
    let (advice_blind, perm_blind) = blind_polys(domain);
    let commitments = round1(&a, &b, &c, advice_blind, &scheme);
    let challenge_generator = ChallengeGenerator::with_digest(&commitments);
    let [beta, lambda] = challenge_generator.generate_challenges();
    let [commit_a, commit_b, commit_c] = commitments;
    let values = advice.map(|e| e.coeffs);

    let permutation_proof = {
        let permutation_proof = circuit.copy_constrains.prove(&values, beta, lambda);
    };
}

fn round1(
    a: &Poly,
    b: &Poly,
    c: &Poly,
    blinding: [Poly; 3],
    scheme: &KzgScheme,
) -> [KzgCommitment; 3] {
    [a, b, c]
        .iter()
        .zip(blinding.iter())
        .map(|(poly, blind)| scheme.commit(*poly + blind))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

//makes blinding polynomials
fn blind_polys(domain: &impl EvaluationDomain<Fr>) -> ([Poly; 3], Poly) {
    let mut rng = rand::thread_rng();
    let scalars = repeat_with(|| Fr::rand(&mut rng))
        .take(9)
        .collect::<Vec<_>>();

    let gates = [
        [scalars[0], scalars[1]],
        [scalars[2], scalars[3]],
        [scalars[4], scalars[5]],
    ];
    let copy = [scalars[6], scalars[7], scalars[8]];

    let gates = gates
        .iter()
        .map(|item| {
            let poly = Poly::from_coefficients_slice(item);
            poly.mul_by_vanishing_poly(domain.clone());
            poly
        })
        .collect::<Vec<_>>();
    let copy = Poly::from_coefficients_slice(&copy);
    (
        gates.try_into().unwrap(),
        copy.mul_by_vanishing_poly(domain.clone()),
    )
}
