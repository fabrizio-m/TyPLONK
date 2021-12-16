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
use plonk::builder::{CircuitBuilder, Variable};

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
    let proof = circuit.prove([3, 4, 6], circuit2);
    assert!(circuit.verify(proof));
}
```