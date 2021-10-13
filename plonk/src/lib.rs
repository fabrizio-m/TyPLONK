use ark_bls12_381::Fr;
use ark_ec::PairingEngine;
use ark_ff::{BigInteger256, FftField, Fp256, UniformRand};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
    Polynomial, UVPolynomial,
};
use commitment::{print_poly, KzgCommitment, KzgScheme};
use permutations::CopyConstrains;
use srs::Srs;
use std::{
    convert::TryInto,
    io::{repeat, Read},
    iter::repeat_with,
};

mod builder;
mod commitment;
mod interpolation;
mod permutations;
mod prof;
mod srs;

pub struct Circuit {
    gate_constrains: GateConstrains,
    copy_constrains: CopyConstrains,
    srs: Srs,
    domain: GeneralEvaluationDomain<Fr>,
}
pub type G1Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G1Affine;
pub type G2Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G2Affine;
pub type Poly = DensePolynomial<Fr>;

struct GateConstrains {
    q_l: Poly,
    q_r: Poly,
    q_0: Poly,
    q_m: Poly,
    q_c: Poly,
}
pub struct Prof {
    a: KzgCommitment,
    b: KzgCommitment,
    c: KzgCommitment,
}

impl Circuit {
    //cosets stolen from dusk
    const CS1: Fr = Self::coset([7, 0, 0, 0]);
    const CS2: Fr = Self::coset([13, 0, 0, 0]);
    const fn coset(parts: [u64; 4]) -> Fr {
        Fr::new(BigInteger256::new(parts))
    }
    pub fn verify(&self, prof: Prof) -> bool {
        true
    }
    //makes a blindign polynomial
    fn blind_polys(n: usize, domain: &impl EvaluationDomain<Fr>) -> ([Poly; 3], Poly) {
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
    fn round1(
        &self,
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
    fn round2(
        &self,
        vars: [&Poly; 3],
        sigmas: [&Poly; 3],
        blinding: Poly,
        challenge: (Fr, Fr),
    ) -> KzgCommitment {
        let vars = Self::iter_evals(self.eval(vars));
        let sigmas = Self::iter_evals(self.eval(sigmas));
        let roots = self.domain.elements();

        let acc = vars
            .into_iter()
            .zip(sigmas.into_iter())
            .zip(roots)
            .map(|((a, b), c)| (a, b, c));

        let (beta, gamma) = challenge;
        let compute = |val, root| val + beta * root + gamma;
        let acc = acc
            .scan(Fr::from(1), |prev, (nom, den, root)| {
                let (n1, n2, n3) = nom;
                let (d1, d2, d3) = den;
                let nominator = compute(n1, root) * compute(n2, root) * compute(n3, root);
                let denominator = compute(d1, root) * compute(d2, root) * compute(d3, root);
                let acc = *prev * (nominator / denominator);
                *prev = acc.clone();
                Some(acc)
            })
            .collect::<Vec<_>>();
        let poly = blinding + Evaluations::from_vec_and_domain(acc, self.domain).interpolate();
        let scheme = KzgScheme::new(&self.srs);
        scheme.commit(poly)
    }
    fn round3(
        &self,
        vars: [&Poly; 3],
        sigmas: [&Poly; 3],
        blinding: Poly,
        challenge: (Fr, Fr),
    ) -> (KzgCommitment, KzgCommitment, KzgCommitment) {
        todo!()
    }
    fn iter_evals(evals: [Evaluations<Fr>; 3]) -> impl Iterator<Item = (Fr, Fr, Fr)> {
        let [a, b, c] = evals;
        let a = a.evals.into_iter();
        let b = b.evals.into_iter();
        let c = c.evals.into_iter();
        a.zip(b).zip(c).map(|((a, b), c)| (a, b, c))
    }
    fn eval(&self, polys: [&Poly; 3]) -> [Evaluations<Fr>; 3] {
        polys
            .iter()
            .map(|ev| ev.evaluate_over_domain_by_ref(self.domain))
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    }
}
#[test]
fn vanish() {
    let poly = Poly::from_coefficients_slice(&[Fr::from(1), Fr::from(2), Fr::from(3)]);
    print_poly(&poly);
    println!("{:#?}", poly);
}
