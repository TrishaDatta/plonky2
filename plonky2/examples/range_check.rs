use anyhow::Result;
use plonky2::field::types::Field;
use plonky2::iop::witness::{PartialWitness, WitnessWrite};
use plonky2::plonk::circuit_builder::CircuitBuilder;
use plonky2::plonk::circuit_data::CircuitConfig;
use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};

/// An example of using Plonky2 to prove that a given value lies in a given range.
fn main() -> Result<()> {
    const D: usize = 2;
    type C = PoseidonGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;

    let config = CircuitConfig::standard_recursion_config();
    let mut builder = CircuitBuilder::<F, D>::new(config);

    let mut values = Vec::new();
    let log_max = 4;

    let len : usize = 300000;

    for i in 0..len {
        let value = builder.add_virtual_target();
        values.push(value);
        builder.range_check(values[i], log_max);
    }

    let mut pw = PartialWitness::new();
    for i in 0..len {
        pw.set_target(values[i], F::from_canonical_usize(3));
    }
    
    
    let data = builder.build::<C>();
    let proof = data.prove(pw)?;

    /*println!(
        "Value {} is less than 2^{}",
        proof.public_inputs[0], log_max,
    );*/

    data.verify(proof)
}
