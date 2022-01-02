use crate::{CompiledCircuit, GateConstrains};
use ark_bls12_381::Fr;
use ark_ff::{One, Zero};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use kgz::{srs::Srs, KzgScheme};
use permutation::{PermutationBuilder, Tag};
use std::{
    collections::HashMap,
    iter::repeat,
    ops::{Add, Mul},
    rc::Rc,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Mutex,
    },
};
#[cfg(test)]
mod test;

#[derive(Debug, Default)]
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
            .find(|size| *size >= rows + 3)
            .unwrap();
        self.gates.resize(size, Gate::Dummy);
    }
    pub fn compile<const I: usize>(circuit: impl Fn([Variable; I])) -> CompiledCircuit<I> {
        //let builder = Rc::new(Mutex::new(Self::new()));
        let context = Context::default();
        let inputs = [(); I].map(|_| Variable::input(&context));
        circuit(inputs);

        let (gates, mut permutation) = context.finish();

        {
            let rows = gates.len();
            let domain = <GeneralEvaluationDomain<Fr>>::new(rows).unwrap();
            let srs = Srs::random(domain.size());
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

#[derive(Hash, Debug, Clone, Copy, PartialEq, Eq)]
pub struct VarId(usize);

#[derive(Default, Debug)]
pub struct InnerContext {
    builder: CircuitBuilder,
    next_var_id: AtomicUsize,
    pending_eq: Vec<(VarId, VarId)>,
    var_map: HashMap<VarId, Tag>,
}

#[derive(Clone, Default, Debug)]
pub struct Context {
    inner: Rc<Mutex<InnerContext>>,
}

impl Context {
    fn new_id(&self) -> VarId {
        let id = self
            .inner
            .lock()
            .unwrap()
            .next_var_id
            .fetch_add(1, Ordering::Relaxed);
        VarId(id)
    }
    fn add_gate(&mut self, gate: Gate) -> usize {
        let builder = &mut self.inner.lock().unwrap().builder;
        builder.add_gate(gate)
    }
    fn add_var(&self, id: VarId, tag: Tag) {
        self.inner.lock().unwrap().var_map.insert(id, tag);
    }
    fn get_var(&self, id: &VarId) -> Option<Tag> {
        self.inner.lock().unwrap().var_map.get(id).cloned()
    }
    fn add_eq(&mut self, left: VarId, right: VarId) {
        println!("new eq: {} = {}", &left.0, &right.0);
        let [a, b] = [left, right].map(|id| self.get_var(&id));
        let context = &mut self.inner.lock().unwrap();
        match a.zip(b) {
            Some((left, right)) => {
                println!("tags: {:?} = {:?}", &left, &right);
                context
                    .builder
                    .permutation
                    .add_constrain(left, right)
                    .unwrap();
            }
            None => {
                context.pending_eq.push((left, right));
            }
        }
    }
    fn finish(mut self) -> (Vec<Gate>, PermutationBuilder<3>) {
        let pending_eq = {
            let mut inner = self.inner.lock().unwrap();
            std::mem::take(&mut inner.pending_eq)
        };
        for (left, right) in pending_eq {
            self.add_eq(left, right);
        }
        let mut inner = Rc::try_unwrap(self.inner).unwrap().into_inner().unwrap();
        assert!(inner.pending_eq.is_empty());
        inner.builder.fill();
        let InnerContext {
            builder: CircuitBuilder {
                gates, permutation, ..
            },
            ..
        } = inner;
        (gates, permutation)
    }
}

#[derive(Clone)]
pub enum Variable {
    Build {
        context: Context,
        id: VarId,
    },
    Compute {
        value: Fr,
        advice_values: Rc<Mutex<[Vec<Fr>; 3]>>,
    },
}
//#[derive(Clone)]
//pub struct Variable2 {
//class: VarClass,
//id: VarId,
//}

impl Variable {
    pub fn assert_eq(&mut self, other: &Variable) {
        match self {
            Variable::Build { context, id } => {
                let left = id;
                match other {
                    Variable::Build { id, .. } => {
                        context.add_eq(*left, *id);
                    }
                    _ => unreachable!(),
                }
            }
            _ => {}
        };
    }
    fn attempt_eq(context: &mut Context, left: (usize, Option<Tag>), right: (usize, Tag)) {}
    fn input(context: &Context) -> Self {
        let id = context.new_id();
        Self::Build {
            id,
            context: context.clone(),
        }
    }

    fn binary_operation(self, right: Variable, operation: GateOperation) -> Variable {
        //let Variable { id, class } = self;
        //let left_id = id;
        //let right_id = right.id;
        match self {
            Variable::Build { mut context, id } => {
                let right_id = match right {
                    Variable::Build { id, .. } => id,
                    _ => {
                        unreachable!()
                    }
                };
                let ids = [(id, 0_usize), (right_id, 1)];

                let j = context.add_gate(operation.build());
                let output_tag = Tag { i: 2, j };
                let id = context.new_id();
                {
                    context.add_var(id, output_tag);
                }
                ids.into_iter()
                    .map(|(id, i)| (id, context.clone().get_var(&id), i))
                    .for_each(|(id, tag, i)| {
                        let mut context = context.clone();
                        match tag {
                            Some(_) => {
                                let new_id = context.new_id();
                                context.add_var(new_id, Tag { i, j });
                                context.add_eq(id, new_id);
                            }
                            None => {
                                context.add_var(id, Tag { i, j });
                            }
                        };
                    });
                Variable::Build {
                    id,
                    context: context.clone(),
                }
            }
            Variable::Compute {
                value,
                advice_values,
            } => {
                let left = value;
                let right = match right {
                    Variable::Compute { value, .. } => value,
                    _ => unreachable!(),
                };
                let value = operation.compute(left, right);
                {
                    let mut advice = advice_values.lock().unwrap();
                    advice
                        .iter_mut()
                        .zip([left, right, value].into_iter())
                        .for_each(|(col, value)| col.push(value));
                }
                Variable::Compute {
                    value,
                    advice_values,
                }
            }
        }
    }
    /*fn is_input(&self) -> bool {
        match self {
            Variable::Build { input, .. } => *input,
            _ => unreachable!(),
        }
    }*/
    fn index(&self) -> usize {
        0
    }
    fn value(&self) -> &Fr {
        match self {
            Variable::Compute { value, .. } => value,
            _ => unreachable!(),
        }
    }
}

enum GateOperation {
    Input,
    Sum,
    Mul,
}

impl GateOperation {
    fn compute(self, a: Fr, b: Fr) -> Fr {
        let d = match self {
            GateOperation::Input => panic!("don't"),
            GateOperation::Sum => a + b,
            GateOperation::Mul => a * b,
        };
        println!("compute result:{}", d);
        d
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
