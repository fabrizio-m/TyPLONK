use crate::CompiledPermutation;
use ark_bls12_381::Fr;
use ark_ff::{One, Zero};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use std::iter::Iterator;

impl<const C: usize> CompiledPermutation<C> {
    pub fn prove(&self, values: &[Vec<Fr>; C], beta: Fr, lambda: Fr) -> Vec<Fr> {
        let cosets = self.cosets;
        let perms = &self.cols;
        assert_eq!(cosets.len(), C);
        let acc = ColIterator::new(values.clone(), perms.clone());
        let acc = acc.scan(Fr::one(), |state, vals| {
            let row_val = vals.into_iter().fold(*state, |state, val| {
                let (cell_val, (tag, value)) = val;
                let numerator = cell_val + beta * tag + lambda;
                let denominator = cell_val + beta * value + lambda;
                let result = numerator / denominator;
                state * result
            });
            Some(row_val)
        });
        let iter_one = std::iter::repeat(Fr::one()).take(1);
        let d = iter_one.clone().chain(acc).chain(iter_one);
        d.collect()
    }
    pub fn verify(
        &self,
        point: usize,
        values: [Fr; C],
        acc_evals: (Fr, Fr),
        beta: Fr,
        lamba: Fr,
    ) -> bool {
        let perms = self.cols.iter().map(|e| e[point].clone());
        let (num, den) = perms
            .zip(values)
            .map(|((label, value), val)| {
                let num = val + beta * label + lamba;
                let den = val + beta * value + lamba;
                (num, den)
            })
            .reduce(|(num1, den1), (num2, den2)| (num1 * num2, den1 * den2))
            .unwrap();
        let rule1 = acc_evals.1 * den - acc_evals.0 * num == Fr::zero();
        rule1
    }
}

struct ColIterator<const C: usize>(
    [std::vec::IntoIter<Fr>; C],
    [std::vec::IntoIter<(Fr, Fr)>; C],
);
impl<const C: usize> ColIterator<C> {
    pub fn new(vals: [Vec<Fr>; C], perms: [Vec<(Fr, Fr)>; C]) -> Self {
        let vals = vals.map(|evals| evals.into_iter());
        let perms = perms.map(|evals| evals.into_iter());
        Self(vals, perms)
    }
}
impl<const C: usize> Iterator for ColIterator<C> {
    type Item = Vec<(Fr, (Fr, Fr))>;

    fn next(&mut self) -> Option<Self::Item> {
        let i1 = self.0.iter_mut();
        let i2 = self.1.iter_mut();
        i1.zip(i2)
            .map(|(val, perm)| {
                let val = val.next()?;
                let perm = perm.next()?;
                Some((val, perm))
            })
            .collect()
    }
}

#[test]
fn domain() {
    use ark_ff::{FpParameters, PrimeField, Zero};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use std::collections::{HashMap, HashSet};

    //println!("T: {}", <<Fr as PrimeField>::Params as FpParameters>::T);
    let domain = <GeneralEvaluationDomain<Fr>>::new(2_usize.pow(16)).unwrap();
    let elems = domain.elements().take(16).for_each(|elem| {
        println!("{}", elem);
    });
    let w = domain.element(1);
    println!("w*w: {}", w * w);
    let set = domain.elements().collect::<HashSet<Fr>>();
    let mut k1 = Fr::from(1);
    while domain.evaluate_vanishing_polynomial(k1).is_zero() {
        k1 += Fr::from(1);
    }
    println!("coset1: {}", k1);
    let coset1 = domain
        .elements()
        .map(|elem| {
            let k = k1 * elem;
            if set.contains(&k) {
                println!("bad one");
            }
            k
        })
        .collect::<HashSet<Fr>>();
    let mut k2 = k1 + Fr::from(1);
    while domain.evaluate_vanishing_polynomial(k1).is_zero() {
        k2 += Fr::from(1);
    }
    let k2 = Fr::from(4);
    println!("coset2: {}", k2);
    let coset2 = domain
        .elements()
        .map(|elem| {
            let k = k2 * elem;
            if set.contains(&k) || coset1.contains(&k) {
                println!("bad one");
            }
            k
        })
        .collect::<HashSet<Fr>>();
}
