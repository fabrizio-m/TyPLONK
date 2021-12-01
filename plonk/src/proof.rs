use crate::{
    builder::Variable,
    utils::{add_to_poly, blind_polys, l0_poly, SlicedPoly},
    CompiledCircuit, GateConstrains, Poly,
};
use ark_bls12_381::Fr;
use ark_ff::{Field, One, Zero};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, Polynomial, UVPolynomial,
};

use challenges::ChallengeGenerator;
use kgz::{KzgCommitment, KzgOpening, KzgScheme};
use permutation::CompiledPermutation;
use std::{
    convert::TryInto,
    fmt::Display,
    ops::{Add, Mul},
    rc::Rc,
    sync::Mutex,
};

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
            //println!("{} {} {}", a[i], b[i], c[i]);
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
    t: [KzgCommitment; 3],
    r: KzgOpening,
}
fn prove<const I: usize>(circuit: &CompiledCircuit<I>, mut advice: [Poly; 3]) -> Proof {
    let scheme = KzgScheme::new(&circuit.srs);
    let domain = &circuit.domain;
    println!("domain size: {}", domain.size());
    let rows = circuit.rows;
    let w = domain.element(1);
    println!("W:{}", w);

    let (advice_blind, perm_blind) = blind_polys(domain);
    let commitments = {
        let [a, b, c] = &advice;
        round1(&a, &b, &c, &scheme)
    };
    let challenge_generator = ChallengeGenerator::with_digest(&commitments);
    let [beta, gamma] = challenge_generator.generate_challenges();
    //let [beta, lambda] = [Fr::from(10), Fr::from(25)];
    println!("beta: {}", beta);
    println!("gamma: {}", gamma);
    let values = advice
        .clone()
        .map(|e| e.evaluate_over_domain(*domain).evals.to_vec());

    let (acc_poly, acc_commitment, acc_poly_w) = {
        let domain = domain;
        let mut evals = circuit.copy_constrains.prove(&values, beta, gamma);
        let acc_shifted = {
            let evals_shifted = &evals[1..];
            Poly::from_coefficients_slice(evals_shifted)
        };
        evals.pop();
        let acc = Evaluations::from_vec_and_domain(evals, domain.clone());
        let acc = acc.interpolate();
        let commitment = scheme.commit(&acc);
        println!("acc degree: {}", acc.degree());
        println!("domain:{}", domain.size());
        let evalsf = acc.clone().evaluate_over_domain(*domain);
        println!();
        println!("accf:");
        for v in evalsf.evals {
            println!("{}", v);
        }
        println!("accw degree: {}", acc_shifted.degree());
        (acc, commitment, acc_shifted)
    };
    let [a, b, c] = advice;
    let mut challenge_generator = ChallengeGenerator::with_digest(&commitments);
    challenge_generator.digest(&acc_commitment);

    let point = challenge_generator.clone().generate_evaluation_point(rows);
    let [alpha] = challenge_generator.generate_challenges();
    println!("point:{}", point);
    let evaluation_point = domain.element(point);
    println!("eval point: {}", evaluation_point);
    let proof = {
        let quotient = quotient_polynomial(
            circuit,
            [&a, &b, &c],
            (&acc_poly, &acc_poly_w),
            [alpha, beta, gamma],
        );
        let quotient_eval = quotient.eval(evaluation_point);
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
        let advice_evals = openings
            .iter()
            .map(|open| open.opening.1)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let [a, b, c] = openings;
        //let evals = acc_poly.evaluate_over_domain(*domain);
        //println!("acc_evals");
        //for e in evals.evals {}
        let z = scheme.open(acc_poly.clone(), evaluation_point);
        let zw = scheme.open(acc_poly.clone(), evaluation_point * w);
        if scheme.verify(&acc_commitment, &z, evaluation_point) {
            println!("z OK");
        }
        if scheme.verify(&acc_commitment, &zw, evaluation_point * w) {
            println!("zw OK");
        }
        println!("z: {}", z.1);
        println!("zw: {}", zw.1);

        let linearisation = linearisation_poly(
            circuit,
            advice_evals,
            [z, zw].map(|open| open.1),
            acc_poly,
            [alpha, beta, gamma],
            evaluation_point,
            &quotient,
        );
        let r = scheme.open(linearisation, evaluation_point);
        let permutation = PermutationProof {
            commitment: acc_commitment,
            z,
            zw,
        };
        let good = r.1 - quotient_eval * domain.evaluate_vanishing_polynomial(evaluation_point);
        assert_eq!(good, Fr::zero());
        let t = quotient.commit(&scheme);
        Proof {
            a,
            b,
            c,
            permutation,
            evaluation_point,
            t,
            r,
        }
    };
    proof
}
fn verify<const I: usize>(circuit: &CompiledCircuit<I>, proof: Proof, scheme: &KzgScheme) -> bool {
    let domain = &circuit.domain;
    let challenges = verify_challenges(&proof, circuit.rows);
    //println!("challenges: ");
    //println!("beta: {}", challenges.0);
    //println!("lambda: {}", challenges.1);
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
        //println!("{}", a);
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
        t,
        r,
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

fn quotient_polynomial<const I: usize>(
    circuit: &CompiledCircuit<I>,
    advice: [&Poly; 3],
    // Z and Zw
    acc: (&Poly, &Poly),
    // challenges alpha, beta, gamma
    challenges: [Fr; 3],
) -> SlicedPoly<3> {
    let domain = &circuit.domain;
    let gates = &circuit.gate_constrains;
    let permutation = &circuit.copy_constrains;
    let GateConstrains {
        q_l,
        q_r,
        q_o,
        q_m,
        q_c,
        ..
    } = gates;
    let [a, b, c] = advice;
    let CompiledPermutation { cols, cosets, .. } = permutation;
    let [alpha, beta, gamma] = challenges;

    let line1 = &(&(q_l.naive_mul(a) + q_r.naive_mul(b)) - &(q_o.naive_mul(c))
        + q_m.naive_mul(a).naive_mul(b))
        + q_c;

    let line2 = [a, b, c]
        .iter()
        .zip(cosets)
        .map(|(advice, coset)| {
            let rhs = DensePolynomial::from_coefficients_vec(vec![*coset * beta, gamma]);
            *advice + &rhs
        })
        .reduce(|one, other| one.naive_mul(&other))
        .unwrap();
    let permutations = cols.clone().map(|col| {
        let eval =
            Evaluations::from_vec_and_domain(col.clone().iter().map(|e| e.1).collect(), *domain);
        eval.interpolate()
    });
    let line3 = [a, b, c]
        .iter()
        .zip(cosets)
        .zip(permutations)
        .map(|((advice, coset), permutation)| {
            let gamma = DensePolynomial::from_coefficients_vec(vec![gamma]);
            let perm =
                permutation.naive_mul(&DensePolynomial::from_coefficients_vec(vec![beta * coset]));
            *advice + &perm + gamma
        })
        .reduce(|one, other| one.naive_mul(&other))
        .unwrap();
    //enforces the first value of acc to be 1
    let line4 = {
        //make lagrange basis l0
        let l0 = l0_poly(*domain);
        let mut acc = acc.0.clone();
        acc.coeffs[0] -= Fr::from(1);
        acc.naive_mul(&l0)
    };
    let one = Fr::one();
    //todo: learn exponenciation
    let combination_element = [one, gamma, gamma, gamma * gamma];
    let constrains = [line1, line2, line3.mul(Fr::from(-1)), line4];
    let target = constrains
        .into_iter()
        .zip(combination_element.into_iter())
        .map(|(constrain, elem)| constrain.mul(elem))
        .reduce(Add::add)
        .unwrap();
    println!("target degree:{}", target.degree());
    {
        println!("eval: {}", target.evaluate(&domain.element(0)));
    };
    //let eval = target.evaluate(&eval_point);
    let target = target.divide_by_vanishing_poly(*domain).unwrap();
    println!("target degree:{}", target.0.degree());
    //println!("target: {:?}", target.1);
    //assert!(target.1.is_zero());
    SlicedPoly::from_poly(target.0, domain.size())
}
fn linearisation_poly<const I: usize>(
    circuit: &CompiledCircuit<I>,
    advice_evals: [Fr; 3],
    acc_evals: [Fr; 2],
    acc: Poly,
    //alpha, beta, gamma
    challenges: [Fr; 3],
    eval_point: Fr,
    t: &SlicedPoly<3>,
) -> Poly {
    let domain = &circuit.domain;
    let gates = &circuit.gate_constrains;
    let permutation = &circuit.copy_constrains;
    let GateConstrains {
        q_l,
        q_r,
        q_o,
        q_m,
        q_c,
        ..
    } = gates;
    let CompiledPermutation { cols, cosets, .. } = permutation;
    let [a, b, c] = advice_evals;
    let [alpha, beta, gamma] = challenges;
    let gate_constrain = q_l.mul(a) + (&(q_r.mul(b)) - &(q_o.mul(c))) + (&q_m.mul(a * b) + q_c);
    let copy_identity = cosets
        .iter()
        .zip(advice_evals.iter())
        .map(|(coset, eval)| *eval + *coset * beta + gamma)
        .reduce(Mul::mul)
        .unwrap();
    let copy_identity = acc.mul(copy_identity);

    let sigma_evals = cols.clone().map(|col| {
        let eval =
            Evaluations::from_vec_and_domain(col.clone().iter().map(|e| e.1).collect(), *domain);
        let poly = eval.interpolate();
        let eval = poly.evaluate(&eval_point);
        (poly, eval)
    });

    let copy_permutation_ab =
        (a + beta * sigma_evals[0].1 + gamma) * (b + beta * sigma_evals[1].1 + gamma);
    let copy_permutation_c = add_to_poly(sigma_evals[2].0.mul(beta), gamma + c);
    let copy_permutation = copy_permutation_c
        .mul(copy_permutation_ab)
        .mul(acc_evals[1]);
    let copy_constrain = &copy_identity - &copy_permutation;
    let l0_eval = l0_poly(*domain);
    let basis_constrain = add_to_poly(acc, Fr::from(-1)).mul(&l0_eval);

    //let combination_element = [one, gamma, gamma * gamma, gamma * gamma * gamma];
    let quotient = t
        .compact(eval_point)
        .mul(domain.evaluate_vanishing_polynomial(eval_point));
    let r = &(gate_constrain + copy_constrain.mul(alpha) + basis_constrain.mul(alpha.square()))
        - &quotient;
    r
}

fn linearisation_commitment<const I: usize>(
    circuit: &CompiledCircuit<I>,
    advice_evals: [Fr; 3],
    acc: KzgCommitment,
    acc_evals: [Fr; 2],
    eval_point: Fr,
    quotient: [KzgCommitment; 3],
    challenges: [Fr; 4],
) -> KzgCommitment {
    let srs = &circuit.srs;
    let scheme = KzgScheme::new(srs);
    let domain = &circuit.domain;
    let cols = domain.size();
    let fixed_commitments = &circuit.gate_constrains.fixed_commitments;
    let permutation = &circuit.copy_constrains;
    let sigma_evals = permutation.sigma_evals(&eval_point);
    let sigma_commitments = permutation.sigma_commitments(&scheme);
    let CompiledPermutation { cosets, .. } = permutation;
    let [alpha, beta, gamma, v] = challenges;

    let line1 = {
        let [a, b, c] = advice_evals;
        let [q_l, q_r, q_o, q_m, q_c] = fixed_commitments;
        q_l * a + q_r * b - q_o * b + q_m * a * b + *q_c
    };
    let line2 = cosets
        .iter()
        .zip(advice_evals.iter())
        .map(|(coset, eval)| *eval + beta * coset + gamma)
        .reduce(Mul::mul)
        .unwrap();
    let l0_eval = l0_poly(*domain).evaluate(&eval_point);
    let line2 = acc * (line2 * alpha + l0_eval * alpha.square());
    let line3 = sigma_evals
        .iter()
        .zip(advice_evals.iter())
        .map(|(sigma_eval, eval)| *eval + beta * sigma_eval + gamma)
        .take(2)
        .reduce(Mul::mul)
        .unwrap();
    let line3 = sigma_commitments[2] * line3 * alpha * beta * acc_evals[1];
    let quotient = <SlicedPoly<3>>::compact_commitment(cols, quotient, eval_point);
    let vanish_eval = domain.evaluate_vanishing_polynomial(eval_point);
    let line4 = quotient * vanish_eval;
    let non_constant = line1 + line2 - line3 - line4;
    let constant_perm = sigma_evals
        .iter()
        .zip(advice_evals.iter())
        .map(|(coset, eval)| *eval + beta * coset + gamma)
        .take(2)
        .reduce(Mul::mul)
        .unwrap();
    let constant_perm = constant_perm * (advice_evals[2] + gamma) * acc_evals[1];
    let constant = -l0_eval * alpha.square() - alpha * constant_perm;
    let scheme = KzgScheme::new(&circuit.srs);

    non_constant - scheme.identity() * constant
}
