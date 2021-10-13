use ark_bls12_381::Fr;
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use std::collections::HashMap;

use crate::Poly;
pub struct CopyConstrains {
    sigma_a: Poly,
    sigma_b: Poly,
    sigma_c: Poly,
}
#[derive(Clone, Copy)]
pub enum Copy {
    Left(usize),
    Right(usize),
    Output(usize),
}

impl From<&HashMap<usize, Vec<Copy>>> for CopyConstrains {
    fn from(copy_constrains: &HashMap<usize, Vec<Copy>>) -> Self {
        let gates = copy_constrains.keys().len();
        let [a, b, c] = Self::permutations(copy_constrains, gates);

        let domain = GeneralEvaluationDomain::new(gates).unwrap();
        let roots = domain.elements().take(gates).collect::<Vec<_>>();

        ///find good ones
        let cosets = (Fr::from(1), Fr::from(2));
        let [sigma_a, sigma_b, sigma_c] = [
            Self::evaluations(a, &roots, cosets),
            Self::evaluations(b, &roots, cosets),
            Self::evaluations(c, &roots, cosets),
        ]
        .map(|evals| {
            let evals = Evaluations::from_vec_and_domain(evals, domain);
            evals.interpolate()
        });
        Self {
            sigma_a,
            sigma_b,
            sigma_c,
        }
    }
}

impl CopyConstrains {
    fn permutations(copy_constrains: &HashMap<usize, Vec<Copy>>, n: usize) -> [Vec<Copy>; 3] {
        let mut a: Vec<Copy> = (0..n).map(|i| Copy::Left(i)).collect();
        let mut b: Vec<Copy> = (0..n).map(|i| Copy::Right(i)).collect();
        let mut c: Vec<Copy> = (0..n).map(|i| Copy::Output(i)).collect();

        for (_gate, constrains) in copy_constrains.iter() {
            for (index, constrain) in constrains.iter().enumerate() {
                let next = match constrains.len() - 1 == index {
                    true => &constrains[0],
                    false => &constrains[index + 1],
                };
                match constrain {
                    Copy::Left(index) => {
                        a[*index] = *next;
                    }
                    Copy::Right(index) => {
                        b[*index] = *next;
                    }
                    Copy::Output(index) => {
                        c[*index] = *next;
                    }
                }
            }
        }
        [a, b, c]
    }
    fn evaluations(permutations: Vec<Copy>, domain: &Vec<Fr>, cosets: (Fr, Fr)) -> Vec<Fr> {
        permutations
            .into_iter()
            .map(|per| match per {
                Copy::Left(index) => domain[index].clone(),
                Copy::Right(index) => cosets.0 * domain[index],
                Copy::Output(index) => cosets.1 * domain[index],
            })
            .collect()
    }
}
