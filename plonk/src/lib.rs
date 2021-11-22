use ark_bls12_381::Fr;
use ark_ec::PairingEngine;
use ark_ff::{UniformRand, Zero};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain,
    Polynomial, UVPolynomial,
};
use kgz::srs::Srs;
use permutation::CompiledPermutation;
use std::{convert::TryInto, iter::repeat_with};

mod builder;
mod proof;

pub type G1Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G1Affine;
pub type G2Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G2Affine;
pub type Poly = DensePolynomial<Fr>;

#[derive(Debug)]
pub struct CompiledCircuit<const INPUTS: usize> {
    gate_constrains: GateConstrains,
    copy_constrains: CompiledPermutation<3>,
    srs: Srs,
    domain: GeneralEvaluationDomain<Fr>,
    pub rows: usize,
}
#[derive(Debug)]
struct GateConstrains {
    q_l: Poly,
    q_r: Poly,
    q_o: Poly,
    q_m: Poly,
    q_c: Poly,
}

impl<const I: usize> CompiledCircuit<I> {
    //makes a blinding polynomial
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
    fn check_row(&self, advice: [Fr; 3], point: Fr) -> bool {
        let constrains = &self.gate_constrains;
        let q_l = constrains.q_l.evaluate(&point);
        let q_r = constrains.q_r.evaluate(&point);
        let q_o = constrains.q_o.evaluate(&point);
        let q_m = constrains.q_m.evaluate(&point);
        let q_c = constrains.q_c.evaluate(&point);
        let [a, b, c] = advice;

        let result = a * q_l + b * q_r - c * q_o + q_m * a * b + q_c;
        //println!("result");
        //println!("{}", result);
        result.is_zero()
    }
}
#[test]
fn vanish() {
    use kgz::print_poly;
    let poly = Poly::from_coefficients_slice(&[Fr::from(1), Fr::from(2), Fr::from(3)]);
    print_poly(&poly);
    println!("{:#?}", poly);
}
