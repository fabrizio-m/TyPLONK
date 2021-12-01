use ark_bls12_381::Fr;
use ark_ff::{BigInteger256, FftField, Fp256, One, UniformRand, Zero};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, MultilinearExtension,
    Polynomial, UVPolynomial,
};
use kgz::{KzgCommitment, KzgScheme, Poly};
use std::{
    collections::{HashMap, HashSet},
    mem::swap,
};

pub mod proving;

#[derive(Hash, Debug, PartialEq, Eq, Clone, Copy)]
pub struct Tag {
    pub i: usize,
    pub j: usize,
}
impl Tag {
    fn to_index(&self, rows: &usize) -> usize {
        self.j + self.i * rows
    }
    fn from_index(index: &usize, rows: &usize) -> Self {
        let i = index / rows;
        let j = index % rows;
        Tag { i, j }
    }
}
#[derive(Default, Debug)]
pub struct PermutationBuilder<const C: usize> {
    //(i,j)
    constrains: HashMap<Tag, Vec<Tag>>,
    rows: usize,
}

impl<const C: usize> PermutationBuilder<C> {
    pub fn add_row(&mut self) {
        self.rows += 1;
    }
    pub fn with_rows(rows: usize) -> Self {
        let mut new = Self::default();
        new.rows = rows;
        new
    }
    ///checks if a tag is valid
    fn check_tag(&self, tag: &Tag) -> bool {
        let Tag { i, j } = tag;
        i <= &C && j < &self.rows
    }
    pub fn add_constrain(&mut self, left: Tag, right: Tag) -> Result<(), ()> {
        if !(self.check_tag(&left) && self.check_tag(&right)) {
            Err(())
        } else {
            let cell = self.constrains.entry(left).or_insert(Vec::with_capacity(8));
            cell.push(right);
            Ok(())
        }
    }
    pub fn add_constrains(&mut self, constrains: Vec<(Tag, Tag)>) {
        for (left, right) in constrains {
            self.add_constrain(left, right).unwrap();
        }
    }
    pub fn build(&mut self) -> Permutation<C> {
        let len = self.rows * C;
        let mut mapping = (0..len).collect::<Vec<_>>();
        let mut aux = (0..len).collect::<Vec<_>>();
        //let mut sizes = (0..len).collect::<Vec<_>>();
        let mut sizes = std::iter::repeat(1).take(len).collect::<Vec<_>>();
        let constrains = std::mem::take(&mut self.constrains);
        //println!("mapping:");
        //print_cycle(&mapping);
        for (left, rights) in constrains.into_iter() {
            let mut left = left.to_index(&self.rows);
            for right in rights {
                let mut right = right.to_index(&self.rows);
                if aux[left] == aux[right] {
                    continue;
                }
                if sizes[aux[left]] < sizes[aux[left]] {
                    swap(&mut left, &mut right);
                }
                sizes[aux[left]] = sizes[aux[left]] + sizes[aux[left]];
                //step 4
                let mut next = right;
                let aux_left = aux[left];
                loop {
                    aux[next] = aux_left;
                    next = mapping[next];
                    if aux[next] == aux_left {
                        break;
                    }
                }
                mapping.swap(left, right);
            }
        }
        //println!("perm:");
        //print_cycle(&mapping);
        Permutation { perm: mapping }
    }
}
#[derive(Debug)]
pub struct Permutation<const C: usize> {
    perm: Vec<usize>,
}

impl<const C: usize> Permutation<C> {
    pub fn compile(self) -> CompiledPermutation<C> {
        assert_eq!(self.perm.len() % C, 0);
        let rows = self.perm.len() / C;
        let cols = self.perm.chunks(rows);
        let cosets = Self::cosets(rows);
        let domain = <GeneralEvaluationDomain<Fr>>::new(rows).unwrap();
        let roots = domain.elements().collect::<Vec<_>>();
        let perm = cols.enumerate().map(|(i, col)| {
            let coefficients = col
                .iter()
                .enumerate()
                .map(|(j, index)| {
                    let tag = Tag::from_index(index, &rows);
                    let value = cosets[tag.i] * roots[tag.j];
                    let tag = cosets[i] * roots[j];
                    (tag, value)
                })
                .collect();
            coefficients
            //let poly = DensePolynomial::from_coefficients_vec(coefficients);
            //poly
        });
        let mut cols: [Vec<(Fr, Fr)>; C] = [0_u8; C].map(|_| Default::default());
        for (i, col) in perm.enumerate() {
            cols[i] = col;
        }
        CompiledPermutation { cols, cosets, rows }
    }
    pub fn print(&self) {
        println!("len: {}", self.perm.len());
        let rows = self.perm.len() / C;
        let perm = &self.perm;
        for j in 0..rows {
            let mut row = vec![];
            for i in 0..C {
                row.push(perm[j + i * rows]);
            }
            println!("{:?}", row);
        }
    }
    fn cosets(gates: usize) -> [Fr; C] {
        let domain = <GeneralEvaluationDomain<Fr>>::new(gates).unwrap();
        let mut cosets = [Fr::zero(); C];

        let mut k = Fr::one();
        for coset in cosets.iter_mut() {
            while domain.evaluate_vanishing_polynomial(k).is_zero() {
                k += Fr::from(1);
            }
            *coset = k;
            k += Fr::from(1);
        }
        cosets
    }
}
#[derive(Debug)]
pub struct CompiledPermutation<const C: usize> {
    //cols: Vec<Vec<(Fr, Fr)>>,
    pub cols: [Vec<(Fr, Fr)>; C],
    pub cosets: [Fr; C],
    rows: usize,
}

impl<const C: usize> CompiledPermutation<C> {
    pub fn sigma_evals(&self, point: &Fr) -> [Fr; C] {
        self.cols
            .iter()
            .map(|col| {
                let poly = col.iter().map(|cell| cell.1).collect();
                let poly = Poly::from_coefficients_vec(poly);
                poly.evaluate(point)
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    }
    pub fn sigma_commitments(&self, scheme: &KzgScheme) -> [KzgCommitment; C] {
        self.cols
            .iter()
            .map(|col| {
                let poly = col.iter().map(|cell| cell.1).collect();
                let poly = Poly::from_coefficients_vec(poly);
                scheme.commit(&poly)
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap()
    }
}

fn print_matrix(matrix: &[usize], cols: usize) {
    println!();
    for row in matrix.chunks(cols) {
        println!("{:?}", row);
    }
}

fn print_cycle(elems: &[usize]) {
    let mut seen = HashSet::new();
    for elem in elems {
        if seen.contains(elem) {
            continue;
        } else {
            seen.insert(*elem);
            let mut cycle = vec![*elem];
            let mut next = elems[*elem];
            loop {
                if seen.contains(&next) {
                    break;
                }
                seen.insert(next);
                cycle.push(next);
                next = elems[next];
            }
            println!("{:?}", cycle);
        }
    }
}
