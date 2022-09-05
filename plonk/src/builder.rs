use crate::{
    description::{CircuitDescription, Var},
    CompiledCircuit, GateConstrains,
};
use ark_bls12_381::Fr;
use ark_ff::{One, Zero};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use kgz::{srs::Srs, KzgScheme};
use permutation::{PermutationBuilder, Tag};
use std::{
    collections::HashMap,
    iter::repeat,
    marker::PhantomData,
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
enum Gate {
    Mul,
    Add,
    //Boolean,
    //Constant(Box<Fr>),
    Dummy,
}

impl CircuitBuilder {
    ///adds a general gate
    fn add_gate(&mut self, gate: Gate) -> usize {
        self.gates.push(gate);
        self.permutation.add_row();
        let index = self.gates.len() - 1;
        index
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

    pub fn compile<const I: usize, C: CircuitDescription<I>>() -> CompiledCircuit<I, C> {
        let circuit = C::run::<BuildVar>;
        let context = Context::default();
        let inputs = [(); I].map(|_| BuildVar::input(&context));
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
                circuit_definition: PhantomData,
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
    fn add_gate(&self, gate: Gate) -> usize {
        let builder = &mut self.inner.lock().unwrap().builder;
        builder.add_gate(gate)
    }
    fn add_var(&self, id: VarId, tag: Tag) {
        self.inner.lock().unwrap().var_map.insert(id, tag);
    }
    fn get_var(&self, id: &VarId) -> Option<Tag> {
        self.inner.lock().unwrap().var_map.get(id).cloned()
    }
    fn add_eq(&self, left: VarId, right: VarId) {
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
    fn finish(self) -> (Vec<Gate>, PermutationBuilder<3>) {
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

impl Variable {
    pub fn assert_eq(&self, other: &Variable) {
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

    fn binary_operation(self, right: Variable, operation: GateOperation) -> Variable {
        match self {
            Variable::Build { context, id } => {
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
                        let context = context.clone();
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
}

enum GateOperation {
    Sum,
    Mul,
}

impl GateOperation {
    fn compute(self, a: Fr, b: Fr) -> Fr {
        let d = match self {
            GateOperation::Sum => a + b,
            GateOperation::Mul => a * b,
        };
        println!("compute result:{}", d);
        d
    }
    fn build(self) -> Gate {
        match self {
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
            Gate::Mul => [Fr::zero(), Fr::zero(), Fr::one(), Fr::one(), Fr::zero()],
            Gate::Add => [Fr::one(), Fr::one(), Fr::one(), Fr::zero(), Fr::zero()],
            Gate::Dummy => [Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero()],
        }
    }
}

#[derive(Clone)]
pub(crate) struct BuildVar {
    context: Context,
    id: VarId,
}
#[derive(Clone)]
pub(crate) struct ComputeVar {
    pub(crate) value: Fr,
    pub(crate) advice_values: Rc<Mutex<[Vec<Fr>; 3]>>,
}

impl BuildVar {
    fn binary_operation(&self, rhs: &Self, operation: GateOperation) -> Self {
        let Self { context, id } = &self;
        let right_id = &rhs.id;
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
                let context = context.clone();
                match tag {
                    Some(_) => {
                        let new_id = context.new_id();
                        context.add_var(new_id, Tag { i, j });
                        context.add_eq(*id, new_id);
                    }
                    None => {
                        context.add_var(*id, Tag { i, j });
                    }
                };
            });
        Self {
            id,
            context: context.clone(),
        }
    }

    fn input(context: &Context) -> Self {
        let id = context.new_id();
        Self {
            id,
            context: context.clone(),
        }
    }
}

impl ComputeVar {
    fn binary_operation(&self, rhs: &Self, operation: GateOperation) -> Self {
        let left = self.value;
        let right = rhs.value;
        let value = operation.compute(left, right);
        {
            let mut advice = self.advice_values.lock().unwrap();
            advice
                .iter_mut()
                .zip([left, right, value].into_iter())
                .for_each(|(col, value)| col.push(value));
        }
        Self {
            value,
            advice_values: self.advice_values.clone(),
        }
    }
}

impl Add for BuildVar {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.binary_operation(&rhs, GateOperation::Sum)
    }
}
impl Add for ComputeVar {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.binary_operation(&rhs, GateOperation::Sum)
    }
}

impl Mul for BuildVar {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.binary_operation(&rhs, GateOperation::Mul)
    }
}
impl Mul for ComputeVar {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.binary_operation(&rhs, GateOperation::Mul)
    }
}
impl Var for BuildVar {
    fn assert_eq(&self, other: &Self) {
        let Self { context, id } = self;
        let left = id;
        context.add_eq(*left, other.id);
    }
}
impl Var for ComputeVar {
    ///would be better to handle the error
    fn assert_eq(&self, _other: &Self) {
        //this would prevent creation of invalid proofs, commented to be able to test invalid proofs
        //assert_eq!(self.value, other.value);
    }
}
