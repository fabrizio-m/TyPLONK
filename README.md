A tiny PLONK implementation.

Todo:
- Maybe add opening batching optimization.
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

    d.assert_eq(&c);
}

fn main() {
    let circuit = CircuitBuilder::compile(circuit);
    let proof = circuit.prove([3, 4, 6], circuit, vec![0]);
    assert!(circuit.verify(proof));
}
```