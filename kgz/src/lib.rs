use crate::{srs::Srs, Poly as MyPoly};
use ark_bls12_381::{Bls12_381, Fr};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::Field;
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use std::{
    fmt::{Debug, Display},
    ops::{Add, Mul, Neg, Sub},
};

pub mod srs;

pub type G1Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G1Affine;
pub type G2Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G2Affine;
pub type Poly = DensePolynomial<Fr>;
pub struct KzgScheme<'a>(&'a Srs);
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KzgCommitment(pub G1Point);

impl KzgCommitment {
    pub fn inner(&self) -> &G1Point {
        &self.0
    }
}
#[derive(Debug, Clone, Copy)]
pub struct KzgOpening(pub G1Point, pub Fr);

impl KzgOpening {
    pub fn eval(self) -> Fr {
        self.1
    }
}

impl<'a> KzgScheme<'a> {
    pub fn new(srs: &'a Srs) -> Self {
        Self(&srs)
    }
    pub fn commit(&self, polynomial: &MyPoly) -> KzgCommitment {
        let commitment = self.evaluate_in_s(polynomial);
        KzgCommitment(commitment)
    }
    fn evaluate_in_s(&self, polynomial: &MyPoly) -> G1Point {
        let srs = self.0.g1_ref();
        assert!(srs.len() > polynomial.degree());
        let poly = polynomial.coeffs.iter();
        let srs = srs.iter();
        let point: G1Point = poly
            .zip(srs)
            .map::<G1Point, _>(|(cof, s)| {
                let d = s.mul(cof.clone());
                d.into()
            })
            .sum();
        point
    }
    pub fn open(&self, mut polynomial: MyPoly, z: impl Into<Fr>) -> KzgOpening {
        let z = z.into();
        let evaluation_at_z = polynomial.evaluate(&z);
        let first = polynomial.coeffs.first_mut().expect("at least 1");
        *first -= evaluation_at_z;
        let root = MyPoly::from_coefficients_slice(&[-(z), 1.into()]);
        let new_poly = &polynomial / &root;
        let opening = self.evaluate_in_s(&new_poly);
        KzgOpening(opening, evaluation_at_z)
    }
    ///verifies the opening P(z) = y
    pub fn verify(
        &self,
        commitment: &KzgCommitment,
        opening: &KzgOpening,
        z: impl Into<Fr> + Debug + Display,
        //y: impl Into<Fr> + Debug,
    ) -> bool {
        let y = opening.1;
        let g1 = self.0.g1_ref();
        let g2s = self.0.g2s_ref();
        let g2 = self.0.g2_ref();
        let a = g2s.clone().into_projective() - (g2.mul(z.into()));
        let b = commitment.0.into_projective() - G1Point::prime_subgroup_generator().mul(y);
        let pairing1 = Bls12_381::pairing(opening.0, a);
        let pairing2 = Bls12_381::pairing(b, g2.clone());
        pairing1 == pairing2
    }
    pub fn identity(&self) -> KzgCommitment {
        let polynomial = Poly::from_coefficients_vec(vec![Fr::from(1)]);
        self.commit(&polynomial)
    }
}
pub fn print_poly(poly: &MyPoly) {
    println!();
    for (i, p) in poly.iter().enumerate() {
        println!("{}.X^{}", p, i);
    }
    println!();
}

#[test]
fn commit() {
    let srs = Srs::from_secret(Fr::from(2), 10);
    let scheme = KzgScheme(&srs);
    let poly = MyPoly::from_coefficients_slice(&[1.into(), 2.into(), 3.into()]);
    let commitment = scheme.commit(&poly);
    let d = Fr::from(1_i32);
    assert_eq!(
        commitment.0.into_projective(),
        G1Point::prime_subgroup_generator().mul(poly.evaluate(&Fr::from(2)))
    );
    assert!(poly.evaluate(&d) == 6.into());
    let opening = scheme.open(poly, d);
    assert!(scheme.verify(&commitment, &opening, d));
}
impl Add for KzgCommitment {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let commitment = self.0 + rhs.0;
        Self(commitment.into())
    }
}
impl Add for KzgOpening {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let eval = self.1 + rhs.1;
        let witness = self.0 + self.0;
        Self(witness, eval)
    }
}
impl Sub for KzgCommitment {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::add(self, -rhs)
    }
}
impl Mul<Fr> for KzgCommitment {
    type Output = Self;

    fn mul(self, rhs: Fr) -> Self::Output {
        let element = self.0.mul(rhs);
        Self(element.into())
    }
}

impl Mul<Fr> for &KzgCommitment {
    type Output = KzgCommitment;

    fn mul(self, rhs: Fr) -> Self::Output {
        let element = self.0.mul(rhs);
        KzgCommitment(element.into())
    }
}
impl Neg for KzgCommitment {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let point = self.0;
        Self(-point)
    }
}
