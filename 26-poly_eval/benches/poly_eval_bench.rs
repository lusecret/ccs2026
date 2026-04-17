#![feature(int_log)]
use criterion::*;
use poly_eval::groups::group::CompressedGroup;
use poly_eval::scalar::Scalar;
use poly_eval::polycommit::pcs_eval::PCSEval;
use poly_eval::groups::group::GroupElement;
use poly_eval::polycommit::pcs_eval::PCSGens;
use poly_eval::groups::group::VartimeMultiscalarMul;
use merlin::Transcript;
use rand::rngs::OsRng;
use rand_core::RngCore;
use std::{thread, time};


fn poly_eval_prove_benchmark(c: &mut Criterion) {
  for &s in [20].iter() { //valid values are 10, 12, 14, 16, 18, 20
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("commits_verify_benchmark");
    group.plot_config(plot_config);

    let num_degrees = (2_usize).pow(s as u32);

    let mut prover_transcript = Transcript::new(b"unit test");
    let polys = PCSGens::new(num_degrees, b"gens_polyeval_verify");
    let mut csprng: OsRng = OsRng;
    let (B, B_RR_vec_c) = PCSEval::prover_setup(&polys.gens.G);

    let B_RR_vec = B_RR_vec_c.iter().map(|C| C.compress()).collect();

    let z_1 = Scalar::random(&mut csprng);

    let a_vec = vec![Scalar::random(&mut csprng); polys.gens.G.len()];
    let z_vec = PCSEval::compute_pow_z(z_1, polys.gens.G.len()); 

    let C = GroupElement::vartime_multiscalar_mul(&a_vec, &polys.gens.G);

    let mut y = Scalar::zero();
    for _index in 0..a_vec.len() {
        y += a_vec[_index] * z_vec[_index];
    }


    let name = format!("ari_circuit_prove_{}", num_degrees);
    group.bench_function(&name, move |b| {
      b.iter(|| {
        let (L_vec, R_vec, B_R_vec,  B_RR_vec_vec, a_opening, z_opening) 
            = PCSEval::prover_prove(&polys.gens.G, &polys.gens.h, &C, &B, &B_RR_vec, z_1, y, &a_vec, &mut prover_transcript);

      });
    });
    group.finish();
  }
}

fn poly_eval_verify_benchmark(c: &mut Criterion) {
  for &s in [20].iter() {  //valid values are 10, 12, 14, 16, 18, 20
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("commits_verify_benchmark");
    group.plot_config(plot_config);

    let num_degrees = (2_usize).pow(s as u32);

    let mut prover_transcript = Transcript::new(b"unit test");
    let polys = PCSGens::new(num_degrees, b"gens_polyeval_verify");
    let mut csprng: OsRng = OsRng;

    let (B, B_RR_vec_c) = PCSEval::prover_setup(&polys.gens.G);

    let B_RR_vec = B_RR_vec_c.iter().map(|C| C.compress()).collect();

    let z_1 = Scalar::random(&mut csprng);

    let a_vec = vec![Scalar::random(&mut csprng); polys.gens.G.len()];
    let z_vec = PCSEval::compute_pow_z(z_1, polys.gens.G.len()); 

    let C = GroupElement::vartime_multiscalar_mul(&a_vec, &polys.gens.G);

    let mut y = Scalar::zero();
    for _index in 0..a_vec.len() {
        y += a_vec[_index] * z_vec[_index];
    }

    let (L_vec, R_vec, B_R_vec,  B_RR_vec_vec, a_opening, z_opening) 
        = PCSEval::prover_prove(&polys.gens.G, &polys.gens.h, &C, &B, &B_RR_vec, z_1, y, &a_vec, &mut prover_transcript);


    let mut verifier_transcript = Transcript::new(b"unit test");

    let (B_delta, B_R, B_RR_vec, k_vec) = PCSEval::verifier_setup(&polys.gens.G);


    let full_B_R_vec = [vec![B_R],  B_R_vec].concat();
    let full_B_RR_vec_vec = [vec![B_RR_vec],  B_RR_vec_vec].concat();
    
    let name = format!("ari_circuit_verify_{}", num_degrees);
    group.bench_function(&name, move |b| {
      b.iter(|| {
        //let mut verifier_transcript = Transcript::new(b"proof");
        let mut verifier_transcript = Transcript::new(b"unit test");
        assert_eq!(
          PCSEval::verifier_verify(&polys.gens.G, &polys.gens.h, &C, &B, &full_B_R_vec, &full_B_RR_vec_vec, z_1, y, 
            &L_vec, &R_vec, &B_delta, &k_vec, &a_opening, &mut verifier_transcript),
          true
        );
      });
    });
    group.finish();
  }
}


fn set_duration() -> Criterion {
  Criterion::default().sample_size(10)
}

criterion_group! {
name = benches_snark;
config = set_duration();
//targets = poly_eval_prove_benchmark
targets = poly_eval_verify_benchmark
}

criterion_main!(benches_snark);
