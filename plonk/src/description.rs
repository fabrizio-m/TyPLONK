use crate::{builder::CircuitBuilder, CompiledCircuit};
use std::ops::{Add, Mul};

pub trait CircuitDescription<const INPUTS: usize>: Sized {
    fn run<V: Var>(inputs: [V; INPUTS]);
    fn build() -> CompiledCircuit<INPUTS, Self> {
        CircuitBuilder::compile::<INPUTS, Self>()
    }
}

pub trait Var
where
    Self: Sized + Add<Output = Self> + Mul<Output = Self> + Clone,
{
    fn assert_eq(&self, other: &Self);
}

pub trait VariableTrait
where
    Self: Sized + Add<Output = Self> + Mul<Output = Self>,
{
}
