A tiny PLONK implementation.

Todo:
- Add support for public inputs.
- Maybe add opening batching optimization.
- Randomize commitments.
- Improve kzg implementation.
- Improve circuit building utilities.
- More examples.
- Benchmarks.
- Improve code and documentation.

```rust 
use plonk::{
    builder::{CircuitBuilder, Variable},
    proof::Proof,
    CompiledCircuit,
};
use ark_bls12_381::Fr;
use ark_poly::EvaluationDomain;

fn circuit(circuit: [Variable; 3]) {
    let [a, b, c] = circuit;
    let a = a.clone() * a;
    let b = b.clone() * b;
    let c = c.clone() * c;
    let mut d = a + b;

    d.equal_to(&c);
}

fn main() {
    let circuit = CircuitBuilder::compile(circuit);
    let proof = circuit.prove([Fr::from(3), Fr::from(4), Fr::from(6)], circuit);
    assert!(circuit.verify(proof));
}
```