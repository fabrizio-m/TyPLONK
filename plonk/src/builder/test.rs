use super::{CircuitBuilder, Variable};
use crate::{proof::Proof, CompiledCircuit};
use ark_bls12_381::Fr;
use ark_poly::EvaluationDomain;
fn circuit1(circuit: [Variable; 5]) {
    let [a, b, c, d, e] = circuit;
    let x = (c + d) + e;
    let mut a = a + b;
    a.equal_to(&x);
}

#[test]
fn circuit1_test() {
    let circuit = CircuitBuilder::compile(circuit1);
    let proof = circuit.prove(
        [
            Fr::from(2),
            Fr::from(7),
            Fr::from(6),
            Fr::from(3),
            Fr::from(4),
        ],
        circuit1,
    );
    println!("proof");
    //println!("{:#?}", proof.a);
    assert!(circuit.verify(proof));
}
