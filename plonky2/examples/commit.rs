use plonky2::field::polynomial::{PolynomialCoeffs, PolynomialValues};
use plonky2::fri::oracle::PolynomialBatch;
use plonky2::plonk::config::{GenericConfig, KeccakGoldilocksConfig};
use plonky2::field::types::Field;
use plonky2::util::timing::TimingTree;
use plonky2::field::fft::fft_root_table;
use plonky2::util::{log2_ceil, log2_strict};
use core::cmp::max;
use std::time::{SystemTime, UNIX_EPOCH};
use proc_status::ProcStatus;

use plonky2::fri::structure::FriInstanceInfo;
use plonky2::fri::structure::FriBatchInfo;
use plonky2::iop::challenger::Challenger;
use plonky2::plonk::plonk_common::PlonkOracle;
use plonky2::fri::structure::FriOracleInfo;
use plonky2::fri::reduction_strategies::FriReductionStrategy;
use plonky2::fri::FriConfig;
use plonky2::fri::structure::FriPolynomialInfo;


// 1 << 17 - 0.33 minutes
// 1 << 20 - 2.65 minutes; 0.24 minutes
// 1 << 24 - 45.240982 minutes

fn print_time_since(start: u128, last: u128, tag: &str) -> u128 {
    let now = SystemTime::now();
    let now_epoc = now
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards");
    let now = now_epoc.as_millis();
    println!("{:?}; time since start {:?}; time since last check: {:?}", tag, (now - start) as f32 / 60000.0, (now - last) as f32 / 60000.0); 
    return now;
}

fn main() {
	const degree : usize = 1 << 8;
	const D: usize = 2;
	// change this
    type C = KeccakGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;

    let x = F::MULTIPLICATIVE_GROUP_GENERATOR;
    let mut vals = Vec::new();
    for _ in 0..degree {
    	vals.push(x);
    }

    let rate_bits = 2;
    let max_quotient_degree_factor = 4;
    let degree_bits = log2_strict(degree);

    let max_fft_points = 1 << (degree_bits + max(rate_bits, log2_ceil(max_quotient_degree_factor)));

    let fft_root_table = fft_root_table(max_fft_points);

	let poly_vals = PolynomialValues::new(vals.clone()).ifft();
    let poly_vals_2 = PolynomialValues::new(vals).ifft();
    let mut values_vec = Vec::new();
    values_vec.push(poly_vals);
    values_vec.push(poly_vals_2);

    let start = SystemTime::now();
    let start_epoch = start
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards");
    let start = start_epoch.as_millis();
    let mut last = start;

    
    let commit = PolynomialBatch::<F, C, D>::from_coeffs(
            values_vec,
            rate_bits,
            false,
            4,
            &mut TimingTree::default(),
            Some(&fft_root_table),
        );

    let mut challenger = Challenger::<F, <KeccakGoldilocksConfig as GenericConfig<D>>::Hasher>::new();

    let zeta = challenger.get_extension_challenge::<D>();

    let zeta_batch : FriBatchInfo<F, D> = FriBatchInfo {
        point: zeta,
        polynomials: FriPolynomialInfo::from_range(0, 0..1),
    };

    let oracles : Vec<FriOracleInfo> = Vec::new();

    let openings = vec![zeta_batch];
    let instance = FriInstanceInfo {
        oracles: oracles,
        batches: openings,
    };

    let fri_config = FriConfig {
        rate_bits: rate_bits,
        cap_height: 4,
        proof_of_work_bits: 16,
        reduction_strategy: FriReductionStrategy::ConstantArityBits(4, 5),
        num_query_rounds: 28,
    };

    let proof = PolynomialBatch::<F, C, D>::prove_openings(
        &instance,
        &[
            &commit,
        ],
        &mut challenger,
        &fri_config.fri_params(degree_bits, false),
        &mut TimingTree::default(),
    );


    let cap = commit.merkle_tree.cap;

    let merkle_caps = &[
        cap
    ];

   /* let fri_challenges = challenger.fri_challenges(merkle_caps, poly_vals, proof.opening_proof.final_poly, degree_bits, fri_config);

    verify_fri_proof::<F, C, D>(
        &instance,
        &proof.openings.to_fri_openings(),
        &fri_challenges,
        merkle_caps,
        &proof.opening_proof,
        &fri_config.fri_params(degree_bits, false),
    )?;
*/



  //  let mem = proc_status::mem_usage().unwrap();
  //  println!("Mem usage in bytes: current={}, peak={}", mem.current, mem.peak);

    
    
}

