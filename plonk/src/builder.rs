use std::collections::HashMap;

use crate::permutations::Copy;
use ark_bls12_381::Fr;
struct CircuitBuilder {
    gates: Vec<Gate>,
    copy_constrains: HashMap<usize, Vec<Copy>>,
}

struct GeneralGate {
    q_l: Fr,
    q_r: Fr,
    q_o: Fr,
    q_m: Fr,
    q_c: Fr,
}
enum Gate {
    General(Box<GeneralGate>),
    Mul,
    Add,
    Boolean,
    Constant(Box<Fr>),
}
struct Var(Option<usize>);

impl CircuitBuilder {
    ///adds a general gate
    fn add_gate(&mut self, gate: Gate, a: Var, b: Var) -> Var {
        self.gates.push(gate);
        let index = self.gates.len() - 1;
        let copy1 = a.0.map(|copy| Copy::Left(copy));
        let copy2 = b.0.map(|copy| Copy::Right(copy));
        let constrains = vec![copy1, copy2].into_iter().filter_map(|e| e).collect();
        self.copy_constrains.insert(index, constrains);
        Var(Some(index))
    }
    ///adds a multiplication gate
    fn mul(&mut self, a: Var, b: Var) -> Var {
        self.add_gate(Gate::Mul, a, b)
    }
    ///adds an addition gate
    fn add(&mut self, a: Var, b: Var) -> Var {
        self.add_gate(Gate::Add, a, b)
    }
    fn add_private_input(&mut self) -> Var {
        Var(None)
    }
}
