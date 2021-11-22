use crate::{CompiledPermutation, PermutationBuilder, Tag};
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
            let mut num = Fr::one();
            let row_val = vals.into_iter().fold(Fr::one(), |state, val| {
                let (cell_val, (tag, value)) = val;
                //println!("cell_val: {}", cell_val);
                let numerator = cell_val + beta * tag + lambda;
                num *= numerator;
                let denominator = cell_val + beta * value + lambda;
                let result = numerator / denominator;
                state * result
            });
            //println!("row_val: {}", row_val);
            //println!("row_num: {}", num);
            *state *= row_val;
            Some(*state)
        });
        let iter_one = std::iter::repeat(Fr::one()).take(1);
        let acc = iter_one.clone().chain(acc).chain(iter_one).collect();
        println!();
        println!("acc: ");
        for v in &acc {
            println!("{}", v);
        }
        acc
    }
    pub fn verify(
        &self,
        point: usize,
        values: [Fr; C],
        acc_evals: (Fr, Fr),
        beta: Fr,
        lambda: Fr,
    ) -> bool {
        println!("acc_eval1: {}", acc_evals.0);
        println!("acc_eval2: {}", acc_evals.1);
        let perms = self.cols.iter().map(|e| e[point].clone());
        let (num, den) = perms
            .zip(values)
            .map(|((label, value), val)| {
                //println!("cell_val: {}", val);
                let num = val + beta * label + lambda;
                let den = val + beta * value + lambda;
                (num, den)
            })
            .reduce(|(num1, den1), (num2, den2)| (num1 * num2, den1 * den2))
            .unwrap();

        println!("num: {}", &num);
        println!("den: {}", &den);

        let lhs = acc_evals.1 * den;
        println!("lhs: {}", lhs);
        let rhs = acc_evals.0 * num;
        println!("rhs: {}", rhs);
        println!("lhs/rhs: {}", lhs / rhs);
        let rule1 = lhs - rhs;
        rule1.is_zero()
    }
    pub fn print(&self, val: bool) {
        let rows = self.rows;
        let cols = &self.cols;
        for j in 0..rows {
            //let mut row = vec![];
            for i in 0..C {
                if val {
                    print!("{}", cols[i][j].1);
                } else {
                    print!("{}", cols[i][j].0);
                }
            }
            println!("");
        }
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
fn perm() {
    let advice = [[1, 2, 3], [4, 5, 2], [7, 8, 1]];
    let mut perm = <PermutationBuilder<3>>::with_rows(3);
    perm.add_constrain(Tag { i: 0, j: 1 }, Tag { i: 1, j: 2 })
        .unwrap();
    perm.add_constrain(Tag { i: 0, j: 0 }, Tag { i: 2, j: 2 })
        .unwrap();
    let (beta, lambda) = (Fr::from(10), Fr::from(25));
    println!("perm: {:?}", perm);
    let perm = perm.build();
    println!("perm builded: ",);
    perm.print();
    let perm = perm.compile();
    println!("perm compiled: ");
    perm.print(false);
    println!("vals");
    perm.print(true);
    let values = advice
        .map(|col| col.map(|v| Fr::from(v)))
        .map(|col| col.to_vec());
    let proof = perm.prove(&values, beta, lambda);
    println!("proof: ");
    for v in &proof {
        println!("{}", v);
    }
    let values = [2, 5, 8].map(|v| Fr::from(v));
    assert!(perm.verify(1, values, (proof[1], proof[2]), beta, lambda));
    let values = [1, 5, 8].map(|v| Fr::from(v));
    assert!(!perm.verify(1, values, (proof[1], proof[2]), beta, lambda));
}