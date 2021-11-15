use crate::{srs::Srs, Poly as MyPoly};
use ark_bls12_381::{Bls12_381, Fr};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use std::fmt::{Debug, Display};

pub mod srs;

pub type G1Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G1Affine;
pub type G2Point = <ark_bls12_381::Bls12_381 as PairingEngine>::G2Affine;
pub type Poly = DensePolynomial<Fr>;
pub struct KzgScheme<'a>(&'a Srs);
#[derive(Debug)]
pub struct KzgCommitment(G1Point);

impl KzgCommitment {
    pub fn inner(&self) -> &G1Point {
        &self.0
    }
}
pub struct KzgOpening(G1Point);

impl<'a> KzgScheme<'a> {
    pub fn new(srs: &'a Srs) -> Self {
        Self(&srs)
    }
    pub fn commit(&self, polynomial: MyPoly) -> KzgCommitment {
        let commitment = self.evaluate_in_s(polynomial);
        KzgCommitment(commitment)
    }
    fn evaluate_in_s(&self, polynomial: MyPoly) -> G1Point {
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
        let opening = self.evaluate_in_s(new_poly);
        KzgOpening(opening)
    }
    pub fn verify(
        &self,
        commitment: KzgCommitment,
        opening: KzgOpening,
        z: impl Into<Fr> + Debug + Display,
        y: impl Into<Fr> + Debug,
    ) -> bool {
        let g1 = self.0.g1_ref();
        let g2s = self.0.g2s_ref();
        let g2 = self.0.g2_ref();
        let a = g2s.clone().into_projective() - (g2.mul(z.into()));
        let b = commitment.0.into_projective() - G1Point::prime_subgroup_generator().mul(y.into());
        let pairing1 = Bls12_381::pairing(opening.0, a);
        let pairing2 = Bls12_381::pairing(b, g2.clone());
        pairing1 == pairing2
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
    let commitment = scheme.commit(poly.clone());
    let d = Fr::from(1_i32);
    assert_eq!(
        commitment.0.into_projective(),
        G1Point::prime_subgroup_generator().mul(poly.evaluate(&Fr::from(2)))
    );
    assert!(poly.evaluate(&d) == 6.into());
    let opening = scheme.open(poly, d);
    assert!(scheme.verify(commitment, opening, d, 6));
}
