use crate::Poly;
use ark_bls12_381::Fr;
use ark_ff::FftField;
use ark_poly::{domain::GeneralEvaluationDomain, EvaluationDomain, UVPolynomial};

#[test]
fn inter() {
    let u = [1.into(), 2.into(), 3.into()];
    let poly = Poly::from_coefficients_slice(&u);
    let domain = <GeneralEvaluationDomain<Fr>>::new(5).unwrap();
    let domain2 = <GeneralEvaluationDomain<Fr>>::new(3).unwrap();
    for e in domain.elements() {
        println!("{}", e);
    }
    println!("domain2:");
    for e in domain2.elements() {
        println!("{}", e);
    }
    println!();
    let mut rng = rand::thread_rng();
    println!("{}", domain.sample_element_outside_domain(&mut rng));
    println!();
    for e in domain.fft(&u) {
        println!("{}", e);
    }
    println!("{:?}", domain.vanishing_polynomial());
    println!("{:?}", poly);
    println!("{:#?}", poly.evaluate_over_domain(domain).interpolate());
}
