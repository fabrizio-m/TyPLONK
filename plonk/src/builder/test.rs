use super::{CircuitBuilder, Variable};

fn circuit1(circuit: [Variable; 5]) {
    let [a, b, c, d, e] = circuit;
    let mut a = a + b;
    let b = (c + d) + e;
    a.equal_to(&b);
}

#[test]
fn circuit1_test() {
    let circuit = CircuitBuilder::compile(circuit1);
    //println!("{:#?}", circuit);
}
