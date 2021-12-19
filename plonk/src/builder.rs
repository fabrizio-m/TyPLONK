use crate::{CompiledCircuit, GateConstrains};
use ark_bls12_381::Fr;
use ark_ff::{One, Zero};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use kgz::{srs::Srs, KzgScheme};
use permutation::{PermutationBuilder, Tag};
use std::{
    iter::repeat,
    ops::{Add, Mul},
    rc::Rc,
    sync::Mutex,
};

#[cfg(test)]
mod test;

#[derive(Debug)]
pub struct CircuitBuilder {
    gates: Vec<Gate>,
    permutation: PermutationBuilder<3>,
}

#[derive(Debug, Clone)]
struct GeneralGate {
    q_l: Fr,
    q_r: Fr,
    q_o: Fr,
    q_m: Fr,
    q_c: Fr,
}
#[derive(Debug, Clone)]
enum Gate {
    General(Box<GeneralGate>),
    Mul,
    Add,
    Boolean,
    Constant(Box<Fr>),
    Dummy,
}

struct Var(Option<usize>);

impl CircuitBuilder {
    fn new() -> Self {
        Self {
            gates: Default::default(),
            permutation: PermutationBuilder::with_rows(0),
        }
    }

    ///adds a general gate
    fn add_gate(&mut self, gate: Gate) -> usize {
        self.gates.push(gate);
        self.permutation.add_row();
        let index = self.gates.len() - 1;
        index
    }
    ///adds a multiplication gate
    fn mul(&mut self) -> usize {
        self.add_gate(Gate::Mul)
    }
    ///adds an addition gate
    fn add(&mut self) -> usize {
        self.add_gate(Gate::Add)
    }
    fn add_private_input(&mut self) -> Var {
        Var(None)
    }
    fn fill(&mut self) {
        let rows = self.gates.len();
        let size = repeat(())
            .scan(2_usize, |state, _| {
                let old = *state;
                *state = old * 2;
                Some(old)
            })
            .find(|size| size >= &rows)
            .unwrap();
        self.gates.resize(size, Gate::Dummy);
    }
    pub fn compile<const I: usize>(circuit: impl Fn([Variable; I])) -> CompiledCircuit<I> {
        let builder = Rc::new(Mutex::new(Self::new()));
        let context = Context { builder };
        let inputs = [(); I].map(|_| Variable::input(&context));
        circuit(inputs);
        let mut builder = Rc::try_unwrap(context.builder)
            .unwrap()
            .into_inner()
            .unwrap();
        builder.fill();

        {
            let rows = builder.gates.len();
            let domain = <GeneralEvaluationDomain<Fr>>::new(rows).unwrap();
            let srs = Srs::random(domain.size());
            let CircuitBuilder {
                gates,
                mut permutation,
            } = builder;
            let mut polys = [(); 5].map(|_| <Vec<Fr>>::with_capacity(rows));
            gates.into_iter().for_each(|gate| {
                let row = gate.to_row();
                polys
                    .iter_mut()
                    .zip(row.into_iter())
                    .for_each(|(col, value)| col.push(value));
            });
            let permutation = permutation.build(rows);
            permutation.print();
            let permutation = permutation.compile();
            let scheme = KzgScheme::new(&srs);
            let polys = polys.map(|evals| {
                let poly = Evaluations::from_vec_and_domain(evals, domain).interpolate();
                let commitment = scheme.commit(&poly);
                (poly, commitment)
            });
            let commitments = polys
                .iter()
                .map(|(_, commitment)| commitment.clone())
                .collect::<Vec<_>>();
            let polys = polys.map(|(poly, _)| poly);

            let [q_l, q_r, q_o, q_m, q_c] = polys;
            let gate_constrains = GateConstrains {
                q_l,
                q_r,
                q_o,
                q_m,
                q_c,
                fixed_commitments: commitments.try_into().unwrap(),
            };
            CompiledCircuit {
                gate_constrains,
                copy_constrains: permutation,
                srs,
                domain,
                rows,
            }
        }
    }
}
#[derive(Clone)]
pub struct Context {
    builder: Rc<Mutex<CircuitBuilder>>,
}

#[derive(Clone)]
pub enum Variable {
    Build {
        context: Context,
        input: bool,
        gate_index: Option<usize>,
    },
    Compute {
        value: Fr,
        advice_values: Rc<Mutex<[Vec<Fr>; 3]>>,
    },
}

impl Variable {
    fn assert_eq(&mut self, other: &Variable) -> bool {
        match self {
            Variable::Build {
                context,
                input,
                gate_index,
            } => {
                let mut builder = context.builder.lock().unwrap();
                builder
                    .permutation
                    .add_constrain(
                        Tag {
                            i: 2,
                            j: gate_index.unwrap(),
                        },
                        Tag {
                            i: 2,
                            j: other.index(),
                        },
                    )
                    .unwrap();
                true
            }
            Variable::Compute {
                value,
                advice_values: _,
            } => value == other.value(),
        }
    }
    fn input(context: &Context) -> Self {
        Self::Build {
            context: context.clone(),
            input: true,
            gate_index: None,
        }
    }

    fn binary_operation(self, right: Variable, operation: GateOperation) -> Variable {
        match self {
            Variable::Build {
                context,
                input,
                gate_index,
            } => {
                let index = {
                    let mut builder = context.builder.lock().unwrap();
                    let left_index = gate_index;
                    let gate_index = builder.add_gate(operation.build());
                    let permutation = &mut builder.permutation;
                    if !input {
                        permutation
                            .add_constrain(
                                Tag {
                                    i: 2,
                                    j: left_index.unwrap(),
                                },
                                Tag {
                                    i: 0,
                                    j: gate_index,
                                },
                            )
                            .unwrap();
                    }
                    if !right.is_input() {
                        permutation
                            .add_constrain(
                                Tag {
                                    i: 2,
                                    j: right.index(),
                                },
                                Tag {
                                    i: 1,
                                    j: gate_index,
                                },
                            )
                            .unwrap();
                    }
                    gate_index
                };
                Variable::Build {
                    context,
                    gate_index: Some(index),
                    input: false,
                }
            }
            Variable::Compute {
                value,
                advice_values,
            } => {
                let a = value;
                let b = *right.value();
                let c = operation.compute(a, b);
                let value = c;
                {
                    let mut advice = advice_values.lock().unwrap();
                    advice
                        .iter_mut()
                        .zip([a, b, c].into_iter())
                        .for_each(|(col, value)| col.push(value));
                }
                Variable::Compute {
                    value,
                    advice_values,
                }
            }
        }
    }
    fn is_input(&self) -> bool {
        match self {
            Variable::Build { input, .. } => *input,
            _ => unreachable!(),
        }
    }
    fn index(&self) -> usize {
        match self {
            Variable::Build { gate_index, .. } => gate_index.unwrap(),
            Variable::Compute { .. } => unreachable!(),
        }
    }
    fn value(&self) -> &Fr {
        match self {
            Variable::Compute { value, .. } => value,
            _ => unreachable!(),
        }
    }
}

enum VarType {
    Input(usize),
    Sum(Rc<Variable>, Rc<Variable>),
    Mul(Rc<Variable>, Rc<Variable>),
}
enum GateOperation {
    Input,
    Sum,
    Mul,
}

impl GateOperation {
    fn compute(self, a: Fr, b: Fr) -> Fr {
        match self {
            GateOperation::Input => panic!("don't"),
            GateOperation::Sum => a + b,
            GateOperation::Mul => a * b,
        }
    }
    fn build(self) -> Gate {
        match self {
            GateOperation::Input => panic!("don't"),
            GateOperation::Sum => Gate::Add,
            GateOperation::Mul => Gate::Mul,
        }
    }
}

impl Add for Variable {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        self.binary_operation(rhs, GateOperation::Sum)
    }
}
impl Mul for Variable {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        self.binary_operation(rhs, GateOperation::Mul)
    }
}

impl Gate {
    fn to_row(self) -> [Fr; 5] {
        match self {
            Gate::General(_) => todo!(),
            Gate::Mul => [Fr::zero(), Fr::zero(), Fr::one(), Fr::one(), Fr::zero()],
            Gate::Add => [Fr::one(), Fr::one(), Fr::one(), Fr::zero(), Fr::zero()],
            Gate::Boolean => todo!(),
            Gate::Constant(_) => todo!(),
            Gate::Dummy => [Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
        }
    }
}
