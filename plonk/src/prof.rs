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
}
struct PermutationProof {
    commitment: KzgCommitment,
    z: KzgOpening,
    zw: KzgOpening,
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
    let w = domain.element(1);

    let (advice_blind, perm_blind) = blind_polys(domain);
    let commitments = {
        let [a, b, c] = &advice;
        round1(&a, &b, &c, advice_blind, &scheme)
    };
    let challenge_generator = ChallengeGenerator::with_digest(&commitments);
    let [beta, lambda] = challenge_generator.generate_challenges();
    //let [commit_a, commit_b, commit_c] = commitments;
    let values = advice.map(|e| e.coeffs);

    let (acc_poly, acc_commitment) = {
        let acc_coeffs = circuit.copy_constrains.prove(&values, beta, lambda);
        let acc = DensePolynomial::from_coefficients_vec(acc_coeffs);
        let acc = acc + perm_blind;
        let commitment = scheme.commit(&acc);
        (acc, commitment)
    };
    let [a, b, c] = values.map(|e| DensePolynomial::from_coefficients_vec(e));
    let mut challenge_generator = ChallengeGenerator::with_digest(&commitments);
    challenge_generator.digest(&acc_commitment);

    let evaluation_point = domain.element(challenge_generator.generate_evaluation_point(rows));
    let proof = {
        //todo: use array zip
        let mut commitments = commitments.into_iter();
        let openings = [a, b, c].map(|poly| {
            let commitment = commitments.next().unwrap();
            let opening = scheme.open(poly, evaluation_point);
            PolyProof {
                commitment,
                opening,
            }
        });
        let [a, b, c] = openings;
        let z = scheme.open(acc_poly.clone(), evaluation_point);
        let zw = scheme.open(acc_poly, evaluation_point * w);
        let permutation = PermutationProof {
            commitment: acc_commitment,
            z,
            zw,
        };
        Proof {
            a,
            b,
            c,
            permutation,
        }
    };
    proof
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
        .map(|(poly, blind)| scheme.commit(&(*poly + blind)))
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
