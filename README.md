A tiny PLONK implementation.

Todo:
- Maybe add opening batching optimization.
- Improve kzg implementation.
- More examples.
- Benchmarks.
- Improve code and documentation.

```rust 

    extern crate plonk;

    use plonk::description::{CircuitDescription, Var};

struct Circuit;
impl CircuitDescription<3> for Circuit {
    fn run<V: Var>(inputs: [V; 3]) {
        let [a, b, c] = inputs;
        let a = a.clone() * a;
        let b = b.clone() * b;
        let c = c.clone() * c;
        let d = a + b;

        d.assert_eq(&c);
    }
}

fn main() {
    let circuit = Circuit::build();
    let proof = circuit.prove([3, 4, 5], vec![0]);
    assert!(circuit.verify(proof));
}
```
