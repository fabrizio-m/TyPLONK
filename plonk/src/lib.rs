use ark_bls12_381::Fr;
use ark_ec::PairingEngine;
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain};
use description::CircuitDescription;
use kgz::{srs::Srs, KzgCommitment};
use permutation::CompiledPermutation;
use std::marker::PhantomData;

pub mod builder;
pub mod description;
mod proof;
mod utils;

pub type G1Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G1Affine;
pub type G2Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G2Affine;
pub type Poly = DensePolynomial<Fr>;

#[derive(Debug)]
pub struct CompiledCircuit<const INPUTS: usize, DESC: CircuitDescription<INPUTS>> {
    gate_constrains: GateConstrains,
    copy_constrains: CompiledPermutation<3>,
    srs: Srs,
    domain: GeneralEvaluationDomain<Fr>,
    circuit_definition: PhantomData<DESC>,
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
