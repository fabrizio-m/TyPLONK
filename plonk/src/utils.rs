use crate::Poly;
use ark_bls12_381::Fr;
use ark_ff::{One, UniformRand, Zero};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, SparsePolynomial},
    EvaluationDomain, Polynomial, UVPolynomial,
};
use kgz::{KzgCommitment, KzgOpening, KzgScheme};
use std::{convert::TryInto, iter::repeat_with, ops::Mul};

pub fn add_to_poly(mut poly: Poly, number: Fr) -> Poly {
    if poly.coeffs.is_empty() {
        Poly::from_coefficients_vec(vec![number])
    } else {
        poly.coeffs[0] += number;
        poly
    }
}

///with [D] for the max degree for each slice, and [S] as the number of slices
#[derive(Debug)]
pub struct SlicedPoly<const S: usize> {
    slices: [Poly; S],
    //the degree of each slice
    degree: usize,
}

impl<const S: usize> SlicedPoly<S> {
    pub fn from_poly(poly: Poly, degree: usize) -> Self {
        assert!(poly.degree() / S <= degree);
        let coeffs = poly.coeffs;
        let mut slices = [(); S].map(|_| Poly::zero());
        coeffs
            .chunks(degree)
            .map(|coeffs| Poly::from_coefficients_slice(coeffs))
            .enumerate()
            .for_each(|(index, slice)| {
                slices[index] = slice;
            });
        Self { slices, degree }
    }
    pub fn commit(&self, scheme: &KzgScheme) -> [KzgCommitment; S] {
        //waiting for array.each_ref()
        self.slices.clone().map(|slice| scheme.commit(&slice))
    }
    pub fn open(&self, scheme: &KzgScheme, point: Fr) -> [KzgOpening; S] {
        self.slices.clone().map(|slice| scheme.open(slice, point))
    }
    pub fn verify_opening(
        commits: &[KzgCommitment; S],
        openings: [KzgOpening; S],
        scheme: &KzgScheme,
        point: Fr,
        degree: usize,
    ) -> Option<Fr> {
        let valid = commits
            .iter()
            .zip(openings.iter())
            .all(|(commitment, opening)| scheme.verify(commitment, opening, point));
        if valid {
            let eval = openings
                .iter()
                .enumerate()
                .map(|(index, slice)| {
                    let exponent =
                        SparsePolynomial::from_coefficients_slice(&[((degree) * index, Fr::one())]);
                    slice.1 * exponent.evaluate(&point)
                })
                .reduce(|one, other| one + other)
                .unwrap();
            Some(eval)
        } else {
            None
        }
    }
    pub fn eval(&self, point: Fr) -> Fr {
        self.slices
            .iter()
            //.rev()
            .enumerate()
            .map(|(index, slice)| {
                let exponent = SparsePolynomial::from_coefficients_slice(&[(
                    (self.degree) * index,
                    Fr::one(),
                )]);
                slice.evaluate(&point) * exponent.evaluate(&point)
            })
            .reduce(|one, other| one + other)
            .unwrap()
    }
    pub fn compact(&self, point: Fr) -> Poly {
        self.slices
            .iter()
            .enumerate()
            .map(|(index, slice)| {
                let exponent = SparsePolynomial::from_coefficients_slice(&[(
                    (self.degree) * index,
                    Fr::one(),
                )]);
                slice.mul(exponent.evaluate(&point))
            })
            .reduce(|one, other| one + other)
            .unwrap()
    }
    pub fn compact_commitment(
        degree: usize,
        commitments: [KzgCommitment; S],
        point: Fr,
    ) -> KzgCommitment {
        commitments
            .iter()
            .enumerate()
            .map(|(index, commit)| {
                let exponent =
                    SparsePolynomial::from_coefficients_slice(&[((degree) * index, Fr::one())]);
                commit * exponent.evaluate(&point)
            })
            .reduce(|one, other| one + other)
            .unwrap()
    }
}

#[test]
fn slicing() {
    use kgz::srs::Srs;
    let eval_point = Fr::from(4);
    let coeffs = [1, 2, 3, 4, 5, 6, 7, 8, 9].map(|e| Fr::from(e));
    let poly = Poly::from_coefficients_slice(&coeffs);
    let eval = poly.evaluate(&eval_point);
    //println!("unsliced eval: {}", eval);
    let sliced_poly = <SlicedPoly<5>>::from_poly(poly, 2);
    //println!("slices:{}", sliced_poly.slices.len());
    //println!("sliced:{:?}", sliced_poly);
    let eval2 = sliced_poly.eval(eval_point);
    //println!("sliced eval: {}", eval);
    assert_eq!(eval, eval2);
    let srs = Srs::random(16);
    let scheme = KzgScheme::new(&srs);
    let commits = sliced_poly.commit(&scheme);
    let openings = sliced_poly.open(&scheme, eval_point);
    let valid = <SlicedPoly<5>>::verify_opening(&commits, openings, &scheme, eval_point, 2);
    assert_eq!(valid.unwrap(), eval);
}

pub fn l0_poly(domain: impl EvaluationDomain<Fr>) -> Poly {
    let n = domain.size();
    let den = Poly::from_coefficients_vec(vec![Fr::from(-1), Fr::from(1)]);
    let den: DenseOrSparsePolynomial<_> = den.mul(Fr::from(n as i32)).into();
    //let vanish = domain.vanishing_polynomial();
    let vanish = SparsePolynomial::from_coefficients_vec(vec![(n, Fr::one()), (0, Fr::from(-1))]);
    let vanish: DenseOrSparsePolynomial<_> = vanish.into();
    let rhs = vanish.divide_with_q_and_r(&den).unwrap().0;
    rhs
}

#[test]
fn l0() {
    use ark_poly::Radix2EvaluationDomain;
    let domain = Radix2EvaluationDomain::new(2_usize.pow(16)).unwrap();
    let l0 = l0_poly(domain);
    println!("w:{}", domain.element(1));
    println!("eval0:{}", l0.evaluate(&Fr::from(0)));
    let evals = l0.evaluate_over_domain(domain);
    println!();
    //for eval in evals.evals.iter() {
    //println!("eval:{}", eval);
    //}
    assert_eq!(
        evals.evals.into_iter().reduce(std::ops::Add::add).unwrap(),
        Fr::one()
    );
}

#[test]
fn sliced() {
    use kgz::srs::Srs;
    let srs = Srs::random(3);
    let scheme = KzgScheme::new(&srs);
    let point = &Fr::from(4);
    let coeffs = [1, 2, 3, 4, 5, 6, 7, 8].map(|e| Fr::from(e));
    let poly = Poly::from_coefficients_slice(&coeffs);
    let whole_eval = poly.evaluate(point);
    println!("whole eval: {}", whole_eval);
    let sliced = <SlicedPoly<3>>::from_poly(poly, 3);
    let sliced_eval = sliced.eval(*point);
    println!("sliced eval: {}", sliced_eval);
    assert_eq!(whole_eval, sliced_eval);
    let compact_poly = sliced.compact(*point);
    let compact_eval = compact_poly.evaluate(point);
    println!("compact eval: {}", compact_eval);
    assert_eq!(whole_eval, compact_eval);
    let compact_commit = scheme.commit(&compact_poly);
    let sliced_commit = sliced.commit(&scheme);
    let sliced_commit = <SlicedPoly<3>>::compact_commitment(3, sliced_commit, *point);
    assert_eq!(compact_commit, sliced_commit);
}
