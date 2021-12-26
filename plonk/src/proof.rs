use crate::{
    builder::Variable,
    utils::{add_to_poly, l0_poly, SlicedPoly},
    CompiledCircuit, GateConstrains, Poly,
};
use ark_bls12_381::Fr;
use ark_ff::{Field, One, UniformRand, Zero};
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
    pub fn prove(
        &self,
        inputs: [impl Into<Fr>; I],
        circuit: impl Fn([Variable; I]),
        public_inputs: Vec<impl Into<Fr>>,
    ) -> Proof {
        let inputs = inputs.map(Into::into);
        let public_inputs = public_inputs
            .into_iter()
            .map(Into::into)
            .collect::<Vec<_>>();

        let advice: [Vec<Fr>; 3] = Default::default();
        let advice = Rc::new(Mutex::new(advice));
        let inputs = inputs.map(|input| Variable::Compute {
            value: input,
            advice_values: advice.clone(),
        });
        circuit(inputs);
        let advice = Rc::try_unwrap(advice).unwrap().into_inner().unwrap();
        let mut rng = rand::thread_rng();
        let advice = advice
            .map(|mut col| {
                col.resize(self.rows - 3, Fr::zero());
                let r = [(); 3].map(|_| Fr::rand(&mut rng));
                col.extend_from_slice(&r);
                col
            })
            .map(|col| Evaluations::from_vec_and_domain(col, self.domain).interpolate());

        let mut public_inputs = public_inputs.to_vec();
        public_inputs.resize(self.rows, Fr::zero());

        let proof = prove(&self, advice, public_inputs);
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
    pub public_inputs: Vec<Fr>,
}
fn prove<const I: usize>(
    circuit: &CompiledCircuit<I>,
    advice: [Poly; 3],
    public_inputs: Vec<Fr>,
) -> Proof {
    let scheme = KzgScheme::new(&circuit.srs);
    let domain = &circuit.domain;
    let w = domain.element(1);

    let public_inputs_poly =
        Evaluations::from_vec_and_domain(public_inputs.clone(), *domain).interpolate();
    let commitments = {
        let [a, b, c] = &advice;
        round1(&a, &b, &c, &scheme)
    };
    let challenge_generator = ChallengeGenerator::with_digest(&commitments);
    let [beta, gamma] = challenge_generator.generate_challenges();
    let values = advice
        .clone()
        .map(|e| e.evaluate_over_domain(*domain).evals.to_vec());

    let (acc_poly, acc_commitment, acc_poly_w) = {
        let domain = domain;
        let mut evals = circuit.copy_constrains.prove(&values, beta, gamma);
        evals.pop();
        let acc_shifted = {
            let mut evals_shifted = evals.clone();
            evals_shifted.rotate_left(1);
            let evals = Evaluations::from_vec_and_domain(evals_shifted, *domain);
            evals.interpolate()
        };
        let acc = Evaluations::from_vec_and_domain(evals, domain.clone());
        let acc = acc.interpolate();
        let commitment = scheme.commit(&acc);
        (acc, commitment, acc_shifted)
    };
    let [a, b, c] = advice;
    let mut challenge_generator = ChallengeGenerator::with_digest(&commitments);
    challenge_generator.digest(&acc_commitment);

    let [alpha, evaluation_point] = challenge_generator.generate_challenges();
    let proof = {
        let public_eval = public_inputs_poly.evaluate(&evaluation_point);
        let quotient = quotient_polynomial(
            circuit,
            [&a, &b, &c],
            (&acc_poly, &acc_poly_w),
            [alpha, beta, gamma],
            public_inputs_poly,
        );
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
        let z = scheme.open(acc_poly.clone(), evaluation_point);
        let zw = scheme.open(acc_poly.clone(), evaluation_point * w);

        let linearisation = linearisation_poly(
            circuit,
            advice_evals,
            [z, zw].map(|open| open.1),
            acc_poly,
            [alpha, beta, gamma],
            evaluation_point,
            &quotient,
            public_eval,
        );
        let r = scheme.open(linearisation, evaluation_point);
        let permutation = PermutationProof {
            commitment: acc_commitment,
            z,
            zw,
        };
        let t = quotient.commit(&scheme);
        Proof {
            a,
            b,
            c,
            permutation,
            evaluation_point,
            t,
            r,
            public_inputs,
        }
    };
    proof
}
fn verify<const I: usize>(circuit: &CompiledCircuit<I>, proof: Proof, scheme: &KzgScheme) -> bool {
    let domain = &circuit.domain;
    let challenges = verify_challenges(&proof);

    let (alpha, beta, gamma, point) = challenges;
    let mut public_inputs = proof.public_inputs.clone();
    public_inputs.resize(circuit.rows, Fr::zero());
    let public_eval = Evaluations::from_vec_and_domain(public_inputs.clone(), *domain)
        .interpolate()
        .evaluate(&point);
    let acc_commitment = proof.permutation.commitment.clone();
    let eval_point = proof.evaluation_point;
    if eval_point != point {
        return false;
    };
    let quotient = proof.t.clone();
    let r_opening = proof.r;
    let (advice, acc) = match verify_openings(proof, scheme, domain.element(1)) {
        Some(evals) => evals,
        None => {
            println!("openings not ok");
            return false;
        }
    };
    let r = linearisation_commitment(
        circuit,
        advice,
        acc_commitment,
        [acc.0, acc.1],
        eval_point,
        quotient,
        [alpha, beta, gamma],
        public_eval,
    );

    let open_valid = scheme.verify(&r, &r_opening, eval_point);
    open_valid && r_opening.1.is_zero()
}
///generates alpha, beta, gamma and the eval point
fn verify_challenges(proof: &Proof) -> (Fr, Fr, Fr, Fr) {
    let commitments = [&proof.a, &proof.b, &proof.c].map(|proof| proof.commitment.clone());
    let challenge_generator = ChallengeGenerator::with_digest(&commitments);
    let challenge: [Fr; 2] = challenge_generator.generate_challenges();
    let mut challenge_generator = ChallengeGenerator::with_digest(&commitments);
    challenge_generator.digest(&proof.permutation.commitment);
    let [alpha, point] = challenge_generator.generate_challenges();

    (alpha, challenge[0], challenge[1], point)
}
fn verify_openings(proof: Proof, scheme: &KzgScheme, w: Fr) -> Option<([Fr; 3], (Fr, Fr))> {
    let Proof {
        a,
        b,
        c,
        permutation,
        evaluation_point,
        ..
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
    public_inputs: Poly,
) -> SlicedPoly<3> {
    let domain = &circuit.domain;
    let w = domain.element(1);
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
        + q_c
        + public_inputs;
    vanishes(&line1, *domain);

    let line2 = [a, b, c]
        .iter()
        .zip(cosets)
        .map(|(advice, coset)| {
            let rhs = DensePolynomial::from_coefficients_vec(vec![gamma, *coset * beta]);
            *advice + &rhs
        })
        .reduce(|one, other| one.naive_mul(&other))
        .unwrap();
    let line2_eval = line2.evaluate(&w);
    let line2 = line2.naive_mul(acc.0);
    let permutations = cols.clone().map(|col| {
        let eval =
            Evaluations::from_vec_and_domain(col.clone().iter().map(|e| e.1).collect(), *domain);
        eval.interpolate()
    });
    let line3 = [a, b, c]
        .iter()
        .zip(permutations)
        .map(|(advice, permutation)| {
            let gamma = DensePolynomial::from_coefficients_vec(vec![gamma]);
            let perm = permutation.mul(beta);
            *advice + &perm + gamma
        })
        .reduce(|one, other| one.naive_mul(&other))
        .unwrap();
    let line3_eval = line3.evaluate(&w);
    assert_eq!(
        line2_eval * acc.0.evaluate(&w) - line3_eval * acc.0.evaluate(&w.square()),
        Fr::zero()
    );
    let line3 = line3.naive_mul(acc.1);
    let line4 = {
        let l0 = l0_poly(*domain);
        let mut acc = acc.0.clone();
        acc.coeffs[0] -= Fr::from(1);
        acc.naive_mul(&l0)
    };
    vanishes(&line4, *domain);
    let combination_element = [Fr::one(), alpha, alpha, alpha.square()];

    let constrains = [line1, line2, -line3, line4];
    let target = constrains
        .into_iter()
        .zip(combination_element.into_iter())
        .map(|(constrain, elem)| constrain.mul(elem))
        .reduce(Add::add)
        .unwrap();

    //vanishes(&target, *domain);
    let target = target.divide_by_vanishing_poly(*domain).unwrap();
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
    public_eval: Fr,
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
    let line1 = q_l.mul(a) + (&(q_r.mul(b)) - &(q_o.mul(c))) + (&q_m.mul(a * b) + q_c);
    let line1 = add_to_poly(line1, public_eval);

    let line2 = cosets
        .iter()
        .zip(advice_evals.iter())
        .map(|(coset, eval)| *eval + *coset * beta * eval_point + gamma)
        .reduce(Mul::mul)
        .unwrap();
    let line2 = acc.mul(line2);

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
    let line3 = copy_permutation_c
        .mul(copy_permutation_ab)
        .mul(acc_evals[1]);

    let copy_constrain = &line2 - &line3;

    let l0_eval = l0_poly(*domain).evaluate(&eval_point);
    let l0_eval = Poly::from_coefficients_vec(vec![l0_eval]);
    let line4 = add_to_poly(acc, Fr::from(-1)).mul(&l0_eval);

    let line5 = t
        .compact(eval_point)
        .mul(domain.evaluate_vanishing_polynomial(eval_point));

    let r = &(line1 + copy_constrain.mul(alpha) + line4.mul(alpha.square())) - &line5;
    r
}

fn linearisation_commitment<const I: usize>(
    circuit: &CompiledCircuit<I>,
    advice_evals: [Fr; 3],
    acc: KzgCommitment,
    acc_evals: [Fr; 2],
    eval_point: Fr,
    quotient: [KzgCommitment; 3],
    //alpha,beta,gamma
    challenges: [Fr; 3],
    public_eval: Fr,
) -> KzgCommitment {
    let srs = &circuit.srs;
    let scheme = KzgScheme::new(srs);
    let domain = &circuit.domain;
    let cols = domain.size();
    let fixed_commitments = &circuit.gate_constrains.fixed_commitments;
    let permutation = &circuit.copy_constrains;
    let sigma_evals = permutation.sigma_evals(&eval_point, *domain);
    let sigma_commitments = permutation.sigma_commitments(&scheme, *domain);
    let CompiledPermutation { cosets, .. } = permutation;
    let [alpha, beta, gamma] = challenges;

    let line1 = {
        let [a, b, c] = advice_evals;
        let [q_l, q_r, q_o, q_m, q_c] = fixed_commitments;
        q_l * a + q_r * b - q_o * c + q_m * a * b + *q_c
    };

    let line2 = cosets
        .iter()
        .zip(advice_evals.iter())
        .map(|(coset, eval)| *eval + beta * coset * eval_point + gamma)
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
    let line5 = quotient * vanish_eval;

    let constant_perm = sigma_evals
        .iter()
        .zip(advice_evals.iter())
        .map(|(sigma_eval, eval)| *eval + beta * sigma_eval + gamma)
        .take(2)
        .reduce(Mul::mul)
        .unwrap();
    let constant_perm = constant_perm * (advice_evals[2] + gamma) * acc_evals[1];
    let constant = alpha * constant_perm + l0_eval * alpha.square() + public_eval;
    let scheme = KzgScheme::new(&circuit.srs);

    line1 + (line2 - (line3 + scheme.identity() * constant)) - line5
}
pub fn vanishes(poly: &Poly, domain: impl EvaluationDomain<Fr>) {
    let (_, rest) = poly.divide_by_vanishing_poly(domain).unwrap();
    assert!(rest.is_zero());
    //println!("vanishes");
}
