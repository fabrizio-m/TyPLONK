use ark_bls12_381::Fr;
use ark_ec::PairingEngine;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain};
use kgz::{srs::Srs, KzgCommitment};
use permutation::CompiledPermutation;

mod builder;
mod proof;
mod utils;

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
    fixed_commitments: [KzgCommitment; 5],
}

#[test]
fn vanish() {
    use ark_poly::UVPolynomial;
    use kgz::print_poly;
    let poly = Poly::from_coefficients_slice(&[Fr::from(1), Fr::from(2), Fr::from(3)]);
    print_poly(&poly);
    println!("{:#?}", poly);
}
