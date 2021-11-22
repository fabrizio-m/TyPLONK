use crate::{builder::Variable, CompiledCircuit, Poly};
use ark_bls12_381::Fr;
use ark_ff::UniformRand;
use ark_poly::{
    domain,
    univariate::{DensePolynomial, SparsePolynomial},
    EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial, Radix2EvaluationDomain,
    UVPolynomial,
};

use challenges::ChallengeGenerator;
use kgz::{KzgCommitment, KzgOpening, KzgScheme};
use std::{convert::TryInto, fmt::Display, iter::repeat_with, rc::Rc, sync::Mutex};

mod challenges;

impl<const I: usize> CompiledCircuit<I> {
    pub fn prove(&self, inputs: [Fr; I], circuit: impl Fn([Variable; I])) -> Proof {
        let advice: [Vec<Fr>; 3] = Default::default();
        let advice = Rc::new(Mutex::new(advice));
        let inputs = inputs.map(|input| Variable::Compute {
            value: input,
            advice_values: advice.clone(),
        });
        circuit(inputs);
        let mut advice = Rc::try_unwrap(advice).unwrap().into_inner().unwrap();
        for c in advice.iter_mut() {
            c.reverse();
            c.reverse();
        }
        for i in 0..self.rows {
            let [a, b, c] = &advice;
            println!("{} {} {}", a[i], b[i], c[i]);
        }
        //let advice = advice.map(|col| DensePolynomial::from_coefficients_vec(col));
        let advice =
            advice.map(|col| Evaluations::from_vec_and_domain(col, self.domain).interpolate());
        let proof = prove(&self, advice);
        proof
    }

    pub fn verify(&self, proof: Proof) -> bool {
        let scheme = KzgScheme::new(&self.srs);
        verify(&self, proof, &scheme)
    }
}

#[derive(Debug)]
pub struct PolyProof {
    commitment: KzgCommitment,
    opening: KzgOpening,
}
impl Display for PolyProof {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "polyprof:")?;
        writeln!(f, "commitment:{:#?}", self.commitment)?;
        writeln!(f, "point:{}", self.opening.0)?;
        writeln!(f, "val:{}", self.opening.1)?;
        Ok(())
    }
}
#[derive(Debug)]
pub struct PermutationProof {
    commitment: KzgCommitment,
    z: KzgOpening,
    zw: KzgOpening,
}
#[derive(Debug)]
pub struct Proof {
    pub a: PolyProof,
    pub b: PolyProof,
    pub c: PolyProof,
    pub permutation: PermutationProof,
    pub evaluation_point: Fr,
}
fn prove<const I: usize>(circuit: &CompiledCircuit<I>, mut advice: [Poly; 3]) -> Proof {
    let scheme = KzgScheme::new(&circuit.srs);
    let domain = &circuit.domain;
    println!("domain size: {}", domain.size());
    /*{
        let p = &advice[0];
        let ev = p.evaluate_over_domain_by_ref(*domain);
        println!("evaluations: ",);
        for e in ev.evals {
            println!("{}", e);
        }
    }*/
    let rows = circuit.rows;
    let w = domain.element(1);
    println!("W:{}", w);

    let (advice_blind, perm_blind) = blind_polys(domain);
    //advice
    //.iter_mut()
    //.zip(advice_blind.iter())
    //.for_each(|(col, blind)| {});
    let commitments = {
        let [a, b, c] = &advice;
        round1(&a, &b, &c, &scheme)
    };
    let challenge_generator = ChallengeGenerator::with_digest(&commitments);
    let [beta, lambda] = challenge_generator.generate_challenges();
    //let [beta, lambda] = [Fr::from(10), Fr::from(25)];
    println!("beta: {}", beta);
    println!("lambda: {}", lambda);
    //let [commit_a, commit_b, commit_c] = commitments;
    let values = advice
        .clone()
        .map(|e| e.evaluate_over_domain(*domain).evals.to_vec());

    let (acc_poly, acc_commitment) = {
        let domain = domain;
        let evals = circuit.copy_constrains.prove(&values, beta, lambda);
        let acc = Evaluations::from_vec_and_domain(evals, domain.clone());
        let acc = acc.interpolate();
        //let acc = acc + perm_blind;
        let commitment = scheme.commit(&acc);
        (acc, commitment)
    };
    //let [a, b, c] = values.map(|e| DensePolynomial::from_coefficients_vec(e));
    let [a, b, c] = advice;
    let mut challenge_generator = ChallengeGenerator::with_digest(&commitments);
    challenge_generator.digest(&acc_commitment);

    let point = challenge_generator.generate_evaluation_point(rows);
    println!("point:{}", point);
    let evaluation_point = domain.element(point);
    println!("eval point: {}", evaluation_point);
    let proof = {
        //todo: use array zip
        let mut commitments = commitments.into_iter();
        let openings = [a, b, c].map(|poly| {
            let commitment = commitments.next().unwrap();
            let opening = scheme.open(poly, evaluation_point);
            PolyProof {
                commitment,
                opening,
            }
        });
        let [a, b, c] = openings;
        //let evals = acc_poly.evaluate_over_domain(*domain);
        //println!("acc_evals");
        //for e in evals.evals {}
        let z = scheme.open(acc_poly.clone(), evaluation_point);
        let zw = scheme.open(acc_poly, evaluation_point * w);
        if scheme.verify(&acc_commitment, &z, evaluation_point) {
            println!("z OK");
        }
        if scheme.verify(&acc_commitment, &zw, evaluation_point * w) {
            println!("zw OK");
        }
        println!("z: {}", z.1);
        println!("zw: {}", zw.1);
        let permutation = PermutationProof {
            commitment: acc_commitment,
            z,
            zw,
        };
        Proof {
            a,
            b,
            c,
            permutation,
            evaluation_point,
        }
    };
    proof
}
fn verify<const I: usize>(circuit: &CompiledCircuit<I>, proof: Proof, scheme: &KzgScheme) -> bool {
    let domain = &circuit.domain;
    let challenges = verify_challenges(&proof, circuit.rows);
    println!("challenges: ");
    println!("beta: {}", challenges.0);
    println!("lambda: {}", challenges.1);
    let (beta, lambda, point) = challenges;
    let (advice, acc) = match verify_openings(proof, scheme, domain.element(1)) {
        Some(evals) => evals,
        None => {
            println!("openings not ok");
            return false;
        }
    };
    println!();
    println!("advice:");
    for a in advice {
        println!("{}", a);
    }
    if !circuit.check_row(advice, domain.element(point)) {
        return false;
    }
    println!("advice ok");
    circuit
        .copy_constrains
        .verify(point, advice, acc, beta, lambda)
}
fn verify_challenges(proof: &Proof, rows: usize) -> (Fr, Fr, usize) {
    let commitments = [&proof.a, &proof.b, &proof.c].map(|proof| proof.commitment.clone());
    let challenge_generator = ChallengeGenerator::with_digest(&commitments);
    let challenge: [Fr; 2] = challenge_generator.generate_challenges();
    let mut challenge_generator = ChallengeGenerator::with_digest(&commitments);
    challenge_generator.digest(&proof.permutation.commitment);
    let perm = challenge_generator.generate_evaluation_point(rows);

    (challenge[0], challenge[1], perm)
    //(Fr::from(10), Fr::from(25), perm)
}
fn verify_openings(proof: Proof, scheme: &KzgScheme, w: Fr) -> Option<([Fr; 3], (Fr, Fr))> {
    let Proof {
        a,
        b,
        c,
        permutation,
        evaluation_point,
    } = proof;
    let advice = [a, b, c];

    let valid = advice.iter().all(|proof| {
        let PolyProof {
            commitment,
            opening,
        } = proof;
        scheme.verify(commitment, opening, evaluation_point)
    });
    if !valid {
        return None;
    };
    let advice = advice.map(|proof| proof.opening.eval());
    let acc = {
        let PermutationProof { commitment, z, zw } = permutation;
        let zv = scheme.verify(&commitment, &z, evaluation_point);
        let zwv = scheme.verify(&commitment, &zw, evaluation_point * w);
        if zv && zwv {
            (z.eval(), zw.eval())
        } else {
            return None;
        }
    };
    Some((advice, acc))
}

fn round1(a: &Poly, b: &Poly, c: &Poly, scheme: &KzgScheme) -> [KzgCommitment; 3] {
    [a, b, c]
        .iter()
        .map(|poly| scheme.commit(poly))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

//makes blinding polynomials
fn blind_polys(domain: &impl EvaluationDomain<Fr>) -> ([Poly; 3], Poly) {
    let mut rng = rand::thread_rng();
    let scalars = repeat_with(|| Fr::rand(&mut rng))
        .take(9)
        .collect::<Vec<_>>();

    let gates = [
        [scalars[0], scalars[1]],
        [scalars[2], scalars[3]],
        [scalars[4], scalars[5]],
    ];
    let copy = [scalars[6], scalars[7], scalars[8]];

    let gates = gates
        .iter()
        .map(|item| {
            let poly = Poly::from_coefficients_slice(item);
            //poly.mul_by_vanishing_poly(domain.clone());
            poly
        })
        .collect::<Vec<_>>();
    let copy = Poly::from_coefficients_slice(&copy);
    (
        gates.try_into().unwrap(),
        copy.mul_by_vanishing_poly(domain.clone()),
    )
}
