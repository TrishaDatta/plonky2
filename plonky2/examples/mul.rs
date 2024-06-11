use anyhow::Result;
use plonky2::field::types::Field;
use plonky2::iop::witness::{PartialWitness, WitnessWrite};
use plonky2::plonk::circuit_builder::CircuitBuilder;
use plonky2::plonk::circuit_data::CircuitConfig;
use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
use std::time::{SystemTime, UNIX_EPOCH};

use rand::rngs::OsRng;
use rand::Rng;

static PIXELS : usize = 1000000;

fn print_time_since(last: u128, tag: &str) -> u128 {
    let now = SystemTime::now();
    let now_epoc = now
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards");
    let now = now_epoc.as_millis();
    println!("{:?} - time since last check: {:?}", tag, (now - last) as f32 / 60000.0); 
    return now;
}

fn main() -> Result<()> {
    const D: usize = 2;
    type C = PoseidonGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;

    let mut rng = OsRng;

    let mut r_vals = Vec::new();

    for i in 0..PIXELS {
        r_vals.push(rng.gen_range(0..256) as u32);
    }
   
     // Timing setup
    let start = SystemTime::now();
    let start_epoch = start
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards");
    let start = start_epoch.as_millis();
    let mut last = start;

    let config = CircuitConfig::standard_recursion_config();
    let mut builder = CircuitBuilder::<F, D>::new(config);

    let mut pw = PartialWitness::new();

    let mut r_targets = Vec::new();
    let mut prod_targets = Vec::new();

    let c = F::from_canonical_u32(30);

    let prod_init = builder.add_virtual_target();
    builder.register_public_input(prod_init);
    prod_targets.push(prod_init);

    for i in 0..PIXELS {
        let r = builder.add_virtual_target();
        r_targets.push(r);

        let p = builder.add_const(r, F::from_canonical_u32(30));
        let mul = builder.mul(p, prod_targets[i]);
        prod_targets.push(mul);
    }

    builder.register_public_input(prod_targets[PIXELS]);

    let data = builder.build::<C>();
    last = print_time_since(last, "setup done"); 

    pw.set_target(prod_targets[0], F::from_canonical_u32(1));

    for i in 0..PIXELS {
        pw.set_target(r_targets[i], F::from_canonical_u32(r_vals[i]));
    }

    let proof = data.prove(pw)?;
    last = print_time_since(last, "proof done"); 

    println!(
        "{}, {}",
        proof.public_inputs[0], proof.public_inputs[1]
    );

  
    let res = data.verify(proof);
    let res = res.unwrap();

    last = print_time_since(last, "verify done"); 

    Ok(())
}
