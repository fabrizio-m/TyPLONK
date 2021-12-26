A tiny PLONK implementation.

Todo:
- Maybe add opening batching optimization.
- Improve kzg implementation.
- More examples.
- Benchmarks.
- Improve code and documentation.

```rust 
use plonk::builder::{CircuitBuilder, Variable};

fn circuit1(circuit: [Variable; 3]) {
    let [a, b, c] = circuit;
    let a = a.clone() * a;
    let b = b.clone() * b;
    let c = c.clone() * c;
    let mut d = a + b;

    d.assert_eq(&c);
}


fn main() {
    let circuit = CircuitBuilder::compile(circuit1);
    let proof = circuit.prove([3, 4, 5], circuit1, vec![0]);
    assert!(circuit.verify(proof));
    let proof = circuit.prove([3, 4, 6], circuit1, vec![0]);
    assert!(!circuit.verify(proof));
}
```