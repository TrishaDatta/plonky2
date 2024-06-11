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
use std::io::{BufRead, BufReader};
use std::fs::File;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Sub;
use rand::{SeedableRng};
use rand_chacha::ChaCha8Rng;

use plonky2_field::goldilocks_field::GoldilocksField;
use plonky2_field::types::Sample;
use rand::rngs::OsRng;

use plonky2::fri::structure::FriInstanceInfo;
use plonky2::fri::structure::FriBatchInfo;
use plonky2::iop::challenger::Challenger;
use plonky2::plonk::plonk_common::PlonkOracle;
use plonky2::fri::structure::FriOracleInfo;
use plonky2::fri::reduction_strategies::FriReductionStrategy;
use plonky2::fri::FriConfig;


static PIXELS : usize = 14;
static EXPONENT : u32 = 3;
static PIXEL_RANGE : i32 = 2_i32.pow(EXPONENT);
static HASH_LENGTH : usize = 128;

fn print_time_since(start: u128, last: u128, tag: &str) -> u128 {
    let now = SystemTime::now();
    let now_epoc = now
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards");
    let now = now_epoc.as_millis();
    println!("{:?}; time since start {:?}; time since last check: {:?}", tag, (now - start) as f32 / 60000.0, (now - last) as f32 / 60000.0); 
    return now;
}

fn get_filename(prefix: &str) -> String {
    let mut filename = prefix.to_owned();
    filename.push_str("image_");
    filename.push_str(&PIXELS.to_string());
    filename.push_str("_");
    filename.push_str(&EXPONENT.to_string());
    filename.push_str(".txt");
    return filename
}

fn read_photo(prefix: &str) -> BufReader<File> {
    let file = File::open(get_filename(prefix)).expect("Unable to open file");
    return BufReader::new(file);
}


fn main() {
    const degree : usize = 1 << 20;
    const D: usize = 2;
    // change this
    type C = KeccakGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;
    let rate_bits = 2;
    let max_quotient_degree_factor = 4;
    let degree_bits = log2_strict(degree);
    let omega = F::primitive_root_of_unity(degree_bits);

    let max_fft_points = 1 << (degree_bits + max(rate_bits, log2_ceil(max_quotient_degree_factor)));

    let fft_root_table = fft_root_table(max_fft_points);

    let mut vals = Vec::new();
    vals.push(F::ONE);
    for _ in 0..degree - 1 {
        vals.push(F::ZERO);
    }
    vals.push(F::ZERO - F::ONE);
    let vanishing_poly = PolynomialCoeffs::new(vals);


    // w_vals = [0, 1,...,PIXEL_RANGE - 1]
    let mut w_vals = Vec::new();
    for i in 0..PIXEL_RANGE {
        let i_in_fr = GoldilocksField(i as u64);
        w_vals.push(i_in_fr);
    }

    let mut w_vals = Vec::new();
    for i in 0..PIXEL_RANGE {
        let i_in_fr = GoldilocksField(i as u64);
        w_vals.push(i_in_fr);
    }

    for i in 0..degree - (PIXEL_RANGE as usize) {
        w_vals.push(F::ZERO);
    }

    // w[X] = poly(w_vals)
    let w = PolynomialValues::new(w_vals).ifft();

    // v_vals = [pixel_0,...,pixel_{D-1}]
    let mut v_vals = Vec::new();
    // z_vals = [sort(v || w)]
    let mut z_vals = Vec::new();

    // reading in photo pixels...
    let file = read_photo("orig_");
    for line in file.lines() {
        let line = line.expect("Unable to read line");
        let i = line.parse::<i32>().unwrap();

        let v_point = GoldilocksField(i as u64); 
        v_vals.push(v_point);
        z_vals.push(i);
    }

    for i in 0..degree - PIXELS {
        v_vals.push(F::ZERO);
    }

    for i in 0..PIXEL_RANGE {
        z_vals.push(i);
    }

    // pad z_vals so that [z[omega*x] - z[x][1 - (z[omega*x] - z[x])] = 0 still holds true
    let z_vals_length = z_vals.len();
    for _ in 0..degree - z_vals_length {
        z_vals.push(PIXEL_RANGE - 1);
    }
    z_vals.sort();

    let mut z_f_vals = Vec::new();
    for i in 0..z_vals.len() {
        z_f_vals.push(GoldilocksField(z_vals[i] as u64));
    }
    
    // v[X] = poly(v_vals)
    let v = PolynomialValues::new(v_vals).ifft();
    // z[X] = poly(z_vals)
    let z = PolynomialValues::new(z_f_vals.clone()).ifft();
    let mut values_vec_0 = Vec::new();
    values_vec_0.push(w.clone());
    values_vec_0.push(v.clone());
    values_vec_0.push(z.clone());

    let start = SystemTime::now();
    let start_epoch = start
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards");
    let start = start_epoch.as_millis();
    let mut last = start;

    
    let commit0 = PolynomialBatch::<F, C, D>::from_coeffs(
            values_vec_0,
            rate_bits,
            false,
            4,
            &mut TimingTree::default(),
            Some(&fft_root_table),
        );

    // REPLACE THIS FIX THIS FIX 
    let mut rng = OsRng;
    let gamma = F::sample(&mut rng);


    // Permutation argument
    // We want to prove:
    //           product_{i=0}^{D-1}(v_i + gamma) * product_{i=0}^{PIXEL_RANGE-1}(w_i + gamma) = product_{i=0}^{D + PIXEL_RANGE - 1}(z_i + gamma) 
    // where v holds the image pixels, w is the range that the pixel values must lie in [0, PIXEL_RANGE-1],
    // and z is the sorted concatentation of v and w

    let mut values_vec_1 = Vec::new();

    // w_prod_vals = [1, (gamma), [(gamma)(1 + gamma)],...,[(gamma)...(PIXEL_RANGE - 1 + gamma)]]
    let mut w_prod_vals = Vec::new();
    let mut product = F::ONE;
    w_prod_vals.push(product);

    for i in 0..PIXEL_RANGE {
        let i_in_fr = GoldilocksField(i as u64);
        product *= i_in_fr + gamma;
        w_prod_vals.push(product);
    }

    let w_prod_vals_len = w_prod_vals.len();
    for _ in 0..degree - w_prod_vals_len {
        product *= gamma;
        w_prod_vals.push(product);
    }
    
    // w_prod_omega_vals = [(gamma), [(gamma)(1 + gamma)],...,[(gamma)...(PIXEL_RANGE + gamma)], 1]
    let mut w_prod_omega_vals = Vec::new();
    for i in 1..w_prod_vals.len() {
        w_prod_omega_vals.push(w_prod_vals[i]);
    }
    w_prod_omega_vals.push(w_prod_vals[0]);

    let w_prod = PolynomialValues::new(w_prod_vals).ifft();

    let w_prod_omega = PolynomialValues::new(w_prod_omega_vals).ifft();

    let mut n_1_coeffs = Vec::new();
    n_1_coeffs.push(omega.exp_u64((degree - 1) as u64));
    n_1_coeffs.push(F::ZERO - F::ONE);
    
    let n_1 = PolynomialCoeffs::from(n_1_coeffs);
    println!("n_1 eval {:?}", n_1.eval(omega.exp_u64((degree - 1) as u64)));

    let mut gamma_coeffs = Vec::new();
    gamma_coeffs.push(gamma);
    let gamma_poly = PolynomialCoeffs::from(gamma_coeffs);
    
    let (q_w, r_w) = (&(&w_prod_omega - &(&w_prod * &(&gamma_poly + &w))) * &n_1).div_rem(&vanishing_poly);
    assert!(r_w.is_zero());

    // Will commit to w_prod[X], q_w[X]
    values_vec_1.push(w_prod);
    values_vec_1.push(q_w);

    // v_prod_vals = [1, (pixel_0 + gamma), [(pixel_0 + gamma)(pixel_1 + gamma)],...,[(pixel_0 + gamma)...(pixel_{D-1} + gamma)]]
    let mut v_prod_vals = Vec::new();
    let mut product = F::ONE;
    v_prod_vals.push(product);

    // reading in photo pixels...
    let file = read_photo("orig_");
    for line in file.lines() {
        let line = line.expect("Unable to read line");
        let i = line.parse::<i32>().unwrap();

        let v_point = GoldilocksField(i as u64); 

        product *= v_point + gamma;
        v_prod_vals.push(product);
    }

    for _ in 0..degree - PIXELS - 1 {
        product *= gamma;
        v_prod_vals.push(product);
    }

    // v_prod_omega_vals = [(pixel_0 + gamma), [(pixel_0 + gamma)(pixel_1 + gamma)],...,[(pixel_0 + gamma)...(pixel_{D-1} + gamma)], 1]
    let mut v_prod_omega_vals = Vec::new();
    for i in 1..v_prod_vals.len() {
        v_prod_omega_vals.push(v_prod_vals[i]);
    }
    v_prod_omega_vals.push(v_prod_vals[0]);

    // for all i \in [1, D + 1], v_prod[omega^i] = \prod_{j=0}^{i-1}(v_j + gamma)
    let v_prod = PolynomialValues::from(v_prod_vals).ifft();

    // v_prod_omega[X] = v_prod[omega*X]
    let v_prod_omega = PolynomialValues::from(v_prod_omega_vals).ifft(); 

    // q_v[X] = (v_prod[omega * X] - (v_prod[X] * (gamma + v[X]))) * n_1[X] / Z_H[X]
    let (q_v, r_v) = (&(&v_prod_omega - &(&v_prod * &(&gamma_poly + &v))) * &n_1).div_rem(&vanishing_poly);
    assert!(r_v.is_zero());

    // Will commit to v_prod[X], q_v[X]
    values_vec_1.push(v_prod);
    values_vec_1.push(q_v);

    // z_prod_vals = [1, z_vals_0 + gamma, [(z_0 + gamma)(z_vals_1 + gamma)],...,[(z_vals_0 + gamma)...(z_vals_{PIXEL_RANGE + D - 1} + gamma)]]
    let mut z_prod_vals = Vec::new();
    let mut product = F::ONE;
    z_prod_vals.push(product);
    for i in 0..z_f_vals.len() - 1 {
        product *= z_f_vals[i] + gamma;
        z_prod_vals.push(product);
    }

    // Range argument
    // We want to prove for the z constructed above that:
    //      (z[X] - z[omega*X])(1 - (z[X] - z[omega*X]) = 0 mod Z_H[X]

    // z_omega_vals = [z_vals_0 + gamma,...,[(z_vals_0 + gamma)...(z_vals_{PIXEL_RANGE + D - 1} + gamma)], 1]
    let mut z_omega_vals = Vec::new();
    for i in 1..z_vals.len() {
        z_omega_vals.push(z_f_vals[i]);
    }
    z_omega_vals.push(z_f_vals[0]);

    // z_prod_omega_vals = [z_vals_0 + gamma, [(z_vals_0 + gamma)(z_vals_1 + gamma)],...,[(z_vals_0 + gamma)...(z_vals_{PIXEL_RANGE + D - 1} + gamma)], 1]
    let mut z_prod_omega_vals = Vec::new();
    for i in 1..z_prod_vals.len() {
        z_prod_omega_vals.push(z_prod_vals[i]);
    }
    z_prod_omega_vals.push(z_prod_vals[0]);

    // for all i \in [1, PIXEL_RANGE + D], z_prod[omega^i] = \prod_{j=0}^{i-1}(z_j + gamma)
    let z_prod = PolynomialValues::from(z_prod_vals).ifft();

    // z_prod_omega[X] = z_prod[omega*X]
    let z_prod_omega = PolynomialValues::from(z_prod_omega_vals).ifft();
    println!("z_omega prods done");

    // q_z[X] = (z_prod[omega * X] - (z_prod[X] * (gamma + z[X]))) * n_1[X] / Z_H[X]
    let (q_z, r_z) = (&(&z_prod_omega - &(&z_prod * &(&gamma_poly + &z))) * &n_1).div_rem(&vanishing_poly);
    assert!(r_z.is_zero());
    println!("q_z prods done");

    let z_omega = PolynomialValues::from(z_omega_vals).ifft();

    let mut one_coeffs = Vec::new();
    one_coeffs.push(F::ONE);
    
    let one = PolynomialCoeffs::from(one_coeffs);

    // q_range[X] = (z[X] - z[omega*X])(1 - (z[X] - z[omega*X]) * n_1[X] / Z_H[X]
    let (q_range, r_range) = (&(&(&z_omega - &z) * &(&one - &(&z_omega - &z))) * &n_1).div_rem_long_division(&vanishing_poly);

    assert!(r_range.is_zero());
    println!("r_range prods done");
    last = print_time_since(start, last, "polynomial commitment done done"); 

    // Will commit to z_prod[X], q_z[X], q_range[X]
    values_vec_1.push(z_prod);
    values_vec_1.push(q_z);
    values_vec_1.push(q_range);

    let commit1 = PolynomialBatch::<F, C, D>::from_coeffs(
            values_vec_1,
            rate_bits,
            false,
            4,
            &mut TimingTree::default(),
            Some(&fft_root_table),
    );


    // Now we prove knowledge of actual hash value (Section 5.5) 
    // Want to generate a[X] and prove that Equation 11 in Section 5.5 holds for
    // this a[X] and the v[X] generated above

    // Use commitments to generate random coefficients [r_0,...,r_{HASH_LENGTH-1}]
    // for random linear combination of sum checks

    // FIX THIS LATER to be hash of things
    let mut hash_coeffs = Vec::new();
    for _ in 0..HASH_LENGTH {
        hash_coeffs.push(F::sample(&mut rng));
    }

    let mut rng = ChaCha8Rng::seed_from_u64(0);
    println!("rng: {:?}", F::sample(&mut rng));

    // Let A be the public hashing matrix (we will generate it with a PRG)
    // a_vals = [\sum_{i=0}{HASH_LENGTH-1}r_i * A_{i, 0},...,\sum_{i=0}{HASH_LENGTH-1}r_i * A_{i, D - 1}]
    let mut a_vals = Vec::new();

    // h_sum_vals = [0, v_vals_0 * a_vals_0 ,..., \sum_{i=0}^{D - 1} v_vals_0 * a_vals_0]
    let mut h_sum_vals = Vec::new();

    // h_sum_omega_vals = [\sum_{i=0}^{1} v_vals_0 * a_vals_0,...,\sum_{i=0}^{D - 1} v_vals_0 * a_vals_0, v_vals_0 * a_vals_0]
    let mut h_sum_omega_vals = Vec::new();
    h_sum_vals.push(F::ZERO);
    let mut sum = F::ZERO;

    // Re-read in pixels
    let file = read_photo("orig_");
    for line in file.lines() {
        let line = line.expect("Unable to read line");
        let i = line.parse::<i32>().unwrap();

        let v_point = GoldilocksField(i as u64); 

        let mut a_point = F::ZERO; 
        for j in 0..hash_coeffs.len() {
            a_point += F::sample(&mut rng) * hash_coeffs[j];
        }
        a_vals.push(a_point);

        sum += v_point * a_point;
        h_sum_vals.push(sum);
        h_sum_omega_vals.push(sum);
    }

    let a_vals_length = a_vals.len();
    for _ in 0..degree - a_vals_length {
        a_vals.push(F::ZERO);
    }

    for _ in 0..degree - PIXELS - 1 {
        h_sum_vals.push(sum);
        h_sum_omega_vals.push(sum);
    }
    h_sum_omega_vals.push(F::ZERO);


    // for all i \in [0, D - 1], a[omega^i] = \sum_{j=0}{HASH_LENGTH-1}r_j * A_{j, i}
    let a = PolynomialValues::from(a_vals).ifft(); 

    // for all i \in [0, D], h_sum[omega^i] = \sum_{j=0}^{i} v_vals_j * a_vals_j
    let h_sum = PolynomialValues::from(h_sum_vals).ifft(); 

    // h_sum_omega[X] = h_sum[omega*X]
    let h_sum_omega = PolynomialValues::from(h_sum_omega_vals).ifft();

    // q_h_sum[X] = (h_sum[omega*X] - h_sum[X] - (v[X] * a[X]))* n_1[X] / Z_H[X]
    let (q_h_sum, r_h_sum) = (&(&(&h_sum_omega - &h_sum) - &(&v * &a))* &n_1).div_rem(&vanishing_poly);
    assert!(r_h_sum.is_zero());
    println!("q_h_sum prods done");

    // Second set of polynomials we commit to
    let mut values_vec_2 = Vec::new();

    // Will commit to a[X], h_sum[X], q_h_sum[X]
    values_vec_2.push(a);
    values_vec_2.push(h_sum);
    values_vec_2.push(q_h_sum);

    let commit2 = PolynomialBatch::<F, C, D>::from_coeffs(
            values_vec_2,
            rate_bits,
            false,
            4,
            &mut TimingTree::default(),
            Some(&fft_root_table),
        );


    last = print_time_since(start, last, "polynomial commitment done done"); 
    println!("poly {:?}", commit0.polynomials.len());
    println!("poly {:?}", commit1.polynomials.len());
    println!("poly {:?}", commit2.polynomials.len());

    let mut challenger = Challenger::<F, <KeccakGoldilocksConfig as GenericConfig<D>>::Hasher>::new();

    let zeta = challenger.get_extension_challenge::<D>();

    let zeta_batch : FriBatchInfo<F, D> = FriBatchInfo {
        point: zeta,
        // fix this NEXT !!!
        polynomials: Vec::new(),
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

    let openings0 = PolynomialBatch::<F, C, D>::prove_openings(
        &instance,
        &[
            &commit0,
            &commit1,
        ],
        &mut challenger,
        &fri_config.fri_params(degree_bits, false),
        &mut TimingTree::default(),
    );

    let cap = commit0.merkle_tree.cap;


    //let mem = proc_status::mem_usage().unwrap();
    //println!("Mem usage in bytes: current={}, peak={}", mem.current, mem.peak);  
    
}

