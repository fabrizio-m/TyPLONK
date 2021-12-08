use crate::{
    builder::{CircuitBuilder, Variable},
    proof::Proof,
    CompiledCircuit,
};
use ark_bls12_381::Fr;
use ark_poly::EvaluationDomain;
fn circuit1(circuit: [Variable; 5]) {
    let [a, b, c, d, e] = circuit;
    let x = (c + d) + e;
    let mut a = a + b;
    a.equal_to(&x);
}
fn circuit2(circuit: [Variable; 3]) {
    let [a, b, c] = circuit;
    let a = a.clone() * a;
    let b = b.clone() * b;
    let c = c.clone() * c;
    let mut d = a + b;

    //let mut a = a + b;
    d.equal_to(&c);
}

#[test]
fn circuit1_test() {
    let circuit = CircuitBuilder::compile(circuit1);
    let proof = circuit.prove(
        [
            Fr::from(2),
            Fr::from(7),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
        ],
        circuit1,
    );
    println!("proof");
    //println!("{:#?}", proof.a);
    assert!(circuit.verify(proof));
}

#[test]
fn circuit2_test() {
    let circuit = CircuitBuilder::compile(circuit2);
    let proof = circuit.prove([Fr::from(3), Fr::from(4), Fr::from(6)], circuit2);
    println!("proof");
    //println!("{:#?}", proof.a);
    assert!(circuit.verify(proof));
}
