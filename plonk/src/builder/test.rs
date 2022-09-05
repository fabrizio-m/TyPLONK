use crate::description::{CircuitDescription, Var};

struct Circuit1;
impl CircuitDescription<5> for Circuit1 {
    fn run<V: Var>(inputs: [V; 5]) {
        let [a, b, c, d, e] = inputs;
        let x = (c + d) + e;
        let a = a + b;
        a.assert_eq(&x);
    }
}
struct Circuit2;
impl CircuitDescription<3> for Circuit2 {
    fn run<V: Var>(inputs: [V; 3]) {
        let [a, b, c] = inputs;
        let a = a.clone() * a;
        let b = b.clone() * b;
        let c = c.clone() * c;
        let d = a + b;

        d.assert_eq(&c);
    }
}

#[test]
fn circuit2_test() {
    let circuit = Circuit2::build();
    let proof = circuit.prove([3, 4, 5], vec![0]);
    assert!(circuit.verify(proof));
}
#[test]
#[should_panic]
fn circuit2_test_bad_inputs() {
    let circuit = Circuit2::build();
    let proof = circuit.prove([3, 4, 6], vec![0]);
    assert!(circuit.verify(proof));
}

#[test]
fn circuit1_test() {
    let circuit = Circuit1::build();
    let proof = circuit.prove([2, 7, 2, 3, 4], vec![0]);
    assert!(circuit.verify(proof));
}
