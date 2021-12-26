use crate::CompiledPermutation;
use ark_bls12_381::Fr;
use ark_ff::One;
use std::iter::Iterator;

impl<const C: usize> CompiledPermutation<C> {
    pub fn prove(&self, values: &[Vec<Fr>; C], beta: Fr, gamma: Fr) -> Vec<Fr> {
        let cosets = self.cosets;
        let perms = &self.cols;
        assert_eq!(cosets.len(), C);
        let acc = ColIterator::new(values.clone(), perms.clone());
        let acc = acc.scan(Fr::one(), |state, vals| {
            let mut num = Fr::one();
            let mut den = Fr::one();
            let row_val = vals.into_iter().fold(Fr::one(), |state, val| {
                let (cell_val, (tag, value)) = val;
                let numerator = cell_val + beta * tag + gamma;
                num *= numerator;
                let denominator = cell_val + beta * value + gamma;
                den *= denominator;
                let result = numerator / denominator;
                state * result
            });

            *state *= row_val;
            Some(*state)
        });
        let iter_one = std::iter::repeat(Fr::one()).take(1);
        let acc = iter_one.clone().chain(acc).collect();
        acc
    }
    pub fn print(&self, val: bool) {
        let rows = self.rows;
        let cols = &self.cols;
        for j in 0..rows {
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
