//use crate::scalar;

use super::super::groups::group::{
    CompressedGroup, CompressedGroupExt, GroupElement, VartimeMultiscalarMul,
    GROUP_BASEPOINT_COMPRESSED,
  };
use super::super::groups::random::RandomTape;
use super::super::groups::transcript::{ProofTranscript,AppendToTranscript};
use super::super::math::Math;
use super::super::scalar::{Scalar, ScalarFromPrimitives};
use super::super::montgomery;
use super::commit_errors::CommitError;
use super::commitments::{Commitments, MultiCommitGens};
//use curve25519_dalek::traits::VartimeMultiscalarMul;
use curve25519_dalek::ristretto::RistrettoPoint;
//use curve25519_dalek::scalar::Scalar;
//use curve25519_dalek::scalar::Scalar;
//use curve25519_dalek::traits::VartimeMultiscalarMul;
//use curve25519_dalek::traits::VartimeMultiscalarMul;
use merlin::Transcript;
use rand_core::OsRng;
use num_bigint::{BigUint, ToBigUint};
use serde_json::to_vec;

#[derive(Debug, PartialEq)]
pub struct PCSGens {
  pub gens: MultiCommitGens, // n+1 G points
}

impl PCSGens {
    pub fn new(nums: usize, label: &'static [u8]) -> PCSGens {
        let gens = MultiCommitGens::new(nums, label);
        PCSGens { gens }
    }

    pub fn new_from(Gs: Vec<GroupElement>, h: GroupElement) -> PCSGens {
        let gens = MultiCommitGens::from(Gs, h);
        PCSGens { gens }
    }
}

#[derive(Debug)]
pub struct PCSEval {
  pub num_degrees: usize, // the number of degrees in the polynomial
  pub Co: Vec<Scalar>, // coefficients of polynomial
}

impl PCSEval {


    pub fn compute_L_R(gens_L: &Vec<GroupElement>, gens_R: &Vec<GroupElement>, u: &GroupElement, 
        a_vec: &Vec<Scalar>, z_vec: &Vec<Scalar>) -> (GroupElement, GroupElement) {
        

        let n_prime = a_vec.len()/2;

        assert_eq!(n_prime, gens_L.len());
        assert_eq!(n_prime, gens_R.len());

        let mut c_L = Scalar::zero();
        let mut c_R = Scalar::zero();
        for _index in 0..n_prime {
            c_L += a_vec[_index] * z_vec[_index + n_prime];
            //println!("------- c_L += {:?} * {:?} for index {:?}", a_vec[_index], z_vec[_index + n_prime], _index);
            c_R += a_vec[_index + n_prime] * z_vec[_index];    
        }

        //println!("------- c_L = {:?} c_R - {:?}", c_L, c_R);

        let mut L = u * c_L;
        let mut R = u * c_R;

        let L_gs: RistrettoPoint = GroupElement::vartime_multiscalar_mul(&a_vec[..n_prime].to_vec(), gens_L);
        let R_gs: RistrettoPoint = GroupElement::vartime_multiscalar_mul(&a_vec[n_prime..].to_vec(), gens_R);
        L += L_gs;
        R += R_gs;

        (L, R)
        //(L_gs, R_gs)
    }

    pub fn compute_pow_z(x: Scalar, n: usize) -> Vec<Scalar> {
        let mut Xx: Vec<Scalar> = Vec::new();
        Xx.push(x);
        for _ in 1..n {
          let ret = Xx.last().unwrap();
          let result = ret.mul(&x);
          Xx.push(result);
        }
        Xx
    }

    pub fn compute_prime_vec(item_vec: &Vec<Scalar>, x: &Scalar) -> Vec<Scalar> {
        let n_prime = item_vec.len()/2;

        let mut ret: Vec<Scalar> = Vec::new();
        for _index in 0..n_prime {
            let item_prime_i = item_vec[_index] + item_vec[_index + n_prime] * x;
            ret.push(item_prime_i);
        }
        ret
    }

    pub fn compute_g_prime_vec(gens: &Vec<GroupElement>, x: &Scalar) -> Vec<GroupElement> {
        let n_prime = gens.len()/2;

        let mut ret: Vec<GroupElement> = Vec::new();
        for _index in 0..n_prime {
            let g_prime_i = gens[_index] + gens[_index + n_prime] * x;
            ret.push(g_prime_i);
        }
        ret
    }

    pub fn compute_halves(gens: &Vec<GroupElement>) -> Vec<GroupElement> {
      //  let mut alpha_prime_vec: Vec<Scalar>  = Vec::new();

        let log_nums = gens.len().ilog2() as usize;
        let mut B_R: Vec<Vec<GroupElement>> = vec![Vec::new(); log_nums];
        

        for _index in 0..gens.len() {
            //let new_part : Scalar  = (0..log_nums).map (|n| {
            for _bit_index in 0..log_nums {
                let bit: usize = (_index >> _bit_index) & 1;
                //println!("try B_R number {:?}, adds index {:?}, size {:?}", log_nums-_bit_index-1, _index, B_R.len());
                if bit == 1 {
                    B_R[log_nums-_bit_index-1].push(gens[_index]);
                  //  println!("B_R number {:?}, adds index {:?}", log_nums-_bit_index-1, _index);
                } 
            }
          //  println!("----------------------");
        } 

        let mut ret: Vec<GroupElement> = Vec::new();
        for _B_R_index in 0..B_R.len() {
           ret.push(B_R[_B_R_index].iter().sum()); 
        }
        ret
               
    }

    pub fn compute_e_vec(z: Scalar, n:usize) -> Vec<Scalar> {
        let mut e_vec: Vec<Scalar> = Vec::new();
        e_vec.push(z);
        for _index in 1..(n.ilog2()) as usize {
            e_vec.push(e_vec[_index-1] * e_vec[_index-1]);
        }
        e_vec        
    }

    pub fn compute_exponents(alpha_vec: &Vec<Scalar>, k_vec: &Vec<Scalar>) -> Vec<Scalar> {
        let n_prime = alpha_vec.len()/2;
        Self::compute_alpha_prime_vec(&alpha_vec[..n_prime].to_vec(), &k_vec[1..].to_vec())
    }

    pub fn compute_alpha_prime_vec(alpha_vec: &Vec<Scalar>, k_vec: &Vec<Scalar>) -> Vec<Scalar> {
        let log_nums = k_vec.len();
        let mut alpha_prime_vec: Vec<Scalar>  = Vec::new();
        for _index in 0..alpha_vec.len() {
            let new_part : Scalar  = (0..k_vec.len()).map (|n| {
                let bit: usize = (_index >> n) & 1;
                //println!("bit is {:?}, index {:?}", bit as usize, log_nums-n-1);
                if bit == 1 {
                    k_vec[log_nums-n-1]
                } else {
                    Scalar::zero()
                }
            }).sum();
            alpha_prime_vec.push(new_part + alpha_vec[_index]);
            //println!("new_part is {:?}", new_part);    
        }
        alpha_prime_vec
    }

    pub fn prover_setup(gens: &Vec<GroupElement>) -> (GroupElement, Vec<GroupElement>) {
        
        let B: GroupElement = gens.iter().sum();
        let B_R_vec = PCSEval::compute_halves(gens);

        (B, B_R_vec)
    }

    pub fn verifier_setup(gens: &Vec<GroupElement>) -> (GroupElement, CompressedGroup, Vec<CompressedGroup>, Vec<Scalar>) {

        let latter_half_gens = &gens[gens.len()/2..].to_vec();
        
        let B: GroupElement = gens.iter().sum();
        let B_R: GroupElement = latter_half_gens.iter().sum();

        let B_R_vec = PCSEval::compute_halves(gens);
        let B_RR_vec = PCSEval::compute_halves(latter_half_gens).iter().map(|C| {C.compress()}).collect();

        let mut setup_transcript = Transcript::new(b"unit test"); // a hack here to use Transcript's challenge generator function

        let k_vec = setup_transcript.challenge_vector(b"generate challenge k_vec", gens.len().ilog2() as usize);

        let B_delta = B + GroupElement::vartime_multiscalar_mul(&k_vec, &B_R_vec);    
        
        (B_delta, B_R.compress(), B_RR_vec, k_vec)
    }

    pub fn prover_prove(gens: &Vec<GroupElement>, u: &GroupElement, C: &GroupElement, B: &GroupElement, B_RR_vec: &Vec<CompressedGroup>,
            z: Scalar, y: Scalar, a_vec: &Vec<Scalar>, transcript: &mut Transcript)
            -> (Vec<CompressedGroup>, Vec<CompressedGroup>, /*3 Vec<CompressedGroup>,*/ Vec<CompressedGroup>, Vec<Vec<CompressedGroup>>, Scalar, Scalar) {
        
        let length = gens.len();
        let length_log2 = gens.len().ilog2() as usize;

        //let x = transcript.challenge_scalar(b"sub-protocol challenge: x");    
        let x = Scalar::one() + Scalar::one() + Scalar::one() + Scalar::one() + Scalar::one();
     //   println!("<<<<< x1 is {:?}", x);

        let u_prime = u * x;
       // let P: GroupElement = C + u_prime * y;   
     //  println!("<<<<< u_primeis {:?}", u_prime.compress());

        let z_vec = Self::compute_pow_z(z, gens.len()); 
        let e_vec = Self::compute_e_vec(z, length);
       
        assert_eq!(e_vec[0], z);
        assert_eq!(z_vec.len(), length);
     //   println!("+++++++++z_vec {:?}", z_vec);

        
        let mut j = 1;

    //3    let alpha_vec = vec![Scalar::one(); length];
        let beta_vec = vec![Scalar::one(); length];

        let mut L_vec: Vec<GroupElement> = Vec::new();
        let mut R_vec: Vec<GroupElement> = Vec::new();
    //3    let mut B_delta_R_vec: Vec<GroupElement> = Vec::new();
        let mut B_R_vec: Vec<CompressedGroup> = Vec::new();   

        let mut B_RR_vec_vec: Vec<Vec<CompressedGroup>> = Vec::new();
     //   B_RR_vec_vec.push(B_RR_vec.clone()); //5 is this needed? since verifier independently can compute this
        //println!("+++++++++B_R_vec_vec.len() {:?}", B_R_vec_vec.len());
      //  println!("+++++++++B_R_vec_vec[0].len() {:?}", B_RR_vec_vec[0].len());
  //      println!("+++++++++B_R_vec_vec.len() {:?}", B_RR_vec_vec.len());
      //  println!("+++++++++B_R_vec_vec[0] {:?}", B_RR_vec_vec[0]);

        Self::prover_prove_recursion(gens, &u_prime, 
            /*3 &alpha_vec,*/ &e_vec, &z_vec, 
            a_vec, L_vec, R_vec, /*B_delta_R_vec,*/B_R_vec, B_RR_vec_vec, transcript)
    }

    pub fn prover_prove_recursion(gens: &Vec<GroupElement>, u: &GroupElement, 
            /*3alpha_vec: &Vec<Scalar>,*/ e_vec: &Vec<Scalar>, z_vec: &Vec<Scalar>, a_vec: &Vec<Scalar>, 
            mut L_vec: Vec<GroupElement>, mut R_vec: Vec<GroupElement>, /* mut B_delta_R_vec: Vec<GroupElement>,*/mut B_R_vec: Vec<CompressedGroup>, mut B_RR_vec_vec: Vec<Vec<CompressedGroup>>, 
            transcript: &mut Transcript)  -> (Vec<CompressedGroup>, Vec<CompressedGroup>, /*3Vec<CompressedGroup>,*/Vec<CompressedGroup>, Vec<Vec<CompressedGroup>>, Scalar, Scalar) {

        if gens.len() == 1 {
            assert_eq!(a_vec.len(),1);
            assert_eq!(z_vec.len(),1);
            //println!("L_vec is {:?} R_vec is {:?} B_delta_R_vec is {:?}", L_vec.len(), R_vec.len(), B_delta_R_vec.len());
    //        println!("a_vec[0] is {:?} z_vec[0] is {:?}", a_vec[0], z_vec[0]);
    //        println!("gens[0] is {:?} ", gens[0].compress());
            let y_c = a_vec[0] * z_vec[0];
            let P_c = gens[0] * a_vec[0] + u * y_c;
    //        println!("----------P_c is {:?}", P_c.compress());            
            let L_vec_c = L_vec.iter().map(|C| {C.compress()}).collect();
            let R_vec_c = R_vec.iter().map(|C| {C.compress()}).collect();
        //3    let B_delta_R_vec_c = B_delta_R_vec.iter().map(|C| {C.compress()}).collect();
            
     //       println!("****B_RR_vec_vec.len() is {:?} ", B_RR_vec_vec.len());

            return (L_vec_c, R_vec_c, /*3B_delta_R_vec_c,*/ B_R_vec, B_RR_vec_vec.to_vec(), a_vec[0], z_vec[0]);
        }
    
        let n_prime =  gens.len()/2;
        let gens_L = gens[n_prime..].to_vec();
        let gens_R = gens[..n_prime].to_vec();

        
     //   println!("enter prover recursion");
        let (L, R) = Self::compute_L_R(&gens_L, &gens_R, u, a_vec, z_vec);


    //1    L.compress().append_to_transcript(b"PCS point L", transcript);
    //1    R.compress().append_to_transcript(b"PCS point R", transcript);

        L_vec.push(L);  
        R_vec.push(R);

    //3    let k_vec = transcript.challenge_vector(b"generate challenge k_vec", gens.len().ilog2() as usize);
    

    //3    let alpha_vec_prime_full = Self::compute_alpha_prime_vec(alpha_vec, &k_vec);
    //3    let alpha_vec_prime = alpha_vec_prime_full[..n_prime].to_vec();

    //3  let B_detal_R = GroupElement::vartime_multiscalar_mul(&alpha_vec_prime_full[n_prime..].to_vec(), &gens_L);
    //1    B_detal_R.compress().append_to_transcript(b"PCS point B_delta_R", transcript);
     //3   B_delta_R_vec.push(B_detal_R);

        let x = transcript.challenge_scalar(b"generate challenge x");
        //let x = Scalar::one()+ Scalar::one() + Scalar::one();
        //let x = Scalar::one();
        let x_inv = x.invert().unwrap();

//        println!("prover x is {:?} x_inv is {:?} ", x, x_inv);

        let gens_vec_prime = Self::compute_g_prime_vec(gens, &x_inv);
        let a_vec_prime = Self::compute_prime_vec(a_vec, &x);
        //4 let B_R_vec_prime = PCSEval::compute_halves(&gens_vec_prime);
        let B_R_prime: &GroupElement = &gens_vec_prime[gens_vec_prime.len() / 2..].iter().sum();
        B_R_vec.push(B_R_prime.compress());

        let B_RR_vec_prime = &PCSEval::compute_halves(&gens_vec_prime[gens_vec_prime.len() / 2..].to_vec());
//        println!("gens_vec_prime/2 is {:?} ", gens_vec_prime.len()/2);  
        
        let z_vec_prime = Self::compute_prime_vec(z_vec, &x_inv);

//        println!("gens_vec_prime is {:?} ", gens_vec_prime.len());  
        if B_RR_vec_prime.len() > 0 {
            B_RR_vec_vec.push(B_RR_vec_prime.iter().map(|C| C.compress()).collect());  
        }
//        println!("B_RR_vec_vec.len() {:?} ", B_RR_vec_vec.len()); 
        Self::prover_prove_recursion(&gens_vec_prime, u, 
            /*&alpha_vec_prime,*/ &e_vec, &z_vec_prime, 
            &a_vec_prime, L_vec, R_vec, /*3B_delta_R_vec,*/ B_R_vec, B_RR_vec_vec, transcript)

    }


    pub fn verifier_verify(gens: &Vec<GroupElement>, u: &GroupElement, C: &GroupElement, B: &GroupElement, B_R_vec: &Vec<CompressedGroup>, B_RR_vec_vec: &Vec<Vec<CompressedGroup>>,
            z: Scalar, y: Scalar, 
            L_vec_c: &Vec<CompressedGroup>, R_vec_c: &Vec<CompressedGroup>, /*B_delta_R_vec_c: &Vec<CompressedGroup>,*/ 
            B_delta: &GroupElement, k_vec: &Vec<Scalar>, a_opening: &Scalar, transcript: &mut Transcript) 
            -> bool {
    

        let L_vec = L_vec_c.iter().map(|C| C.decompress().unwrap()).collect();
        let R_vec = R_vec_c.iter().map(|C| C.decompress().unwrap()).collect();
     //3   let B_delta_R_vec = B_delta_R_vec_c.iter().map(|C| C.decompress().unwrap()).collect();

        let length = gens.len();
        let length_log2 = gens.len().ilog2() as usize;

        //let x = transcript.challenge_scalar(b"sub-protocol challenge: x");  
        let x = Scalar::one() + Scalar::one() + Scalar::one() + Scalar::one() + Scalar::one();  
   //     println!(">>>>> x1 is {:?}", x);

        let e_vec = Self::compute_e_vec(z, length);
       
     //   println!(">>>>> e_vec is {:?}", e_vec);
        assert_eq!(e_vec[0], z);
        //assert_eq!(e_vec.len(), length_log2-1);
        //assert_eq!(z_vec.len(), length);
        
        let z_l = z;
        let z_r = z * e_vec[e_vec.len()-1];

        let u_prime = u * x;
        let P: GroupElement = C + u_prime *y;
       
        let j = 0usize;
        let mut beta_vec = vec![Scalar::zero(); length_log2];

   //     println!(">>>>> u_primeis {:?}", u_prime.compress());

   let B_R_vec_decompressed: Vec<GroupElement> = B_R_vec.iter().map(|C| C.decompress().unwrap()).collect();

        Self::verifier_verify_recursion(length, &u_prime, &P, B, /*1&B.clone(), */
            z_l, z_r, y, j, &e_vec, &mut beta_vec,
            &L_vec, &R_vec, /*3  &B_delta_R_vec,*/ &B_R_vec_decompressed, B_RR_vec_vec, B_delta, k_vec, a_opening, transcript) 
    }

    pub fn verifier_verify_recursion(n: usize, u_update: &GroupElement, P: &GroupElement, B: &GroupElement, /*2B_delta: &GroupElement,*/ 
            z_l: Scalar, z_r: Scalar, y: Scalar, j: usize, e_vec: &Vec<Scalar>, beta_vec: &mut Vec<Scalar>,
            L_vec: &Vec<GroupElement>, R_vec: &Vec<GroupElement>, /*3 B_delta_R_vec: &Vec<GroupElement>,*/B_R_vec: &Vec<GroupElement>, B_RR_vec_vec: &Vec<Vec<CompressedGroup>>,
            B_delta: &GroupElement, k_vec: &Vec<Scalar>, a_opening: &Scalar, transcript: &mut Transcript) 
            -> bool {
            //    println!("round  j ================================================================ {:?}", j);
        if n == 1 {
           // println!("a_opening is {:?}", a_opening);
 //           println!("B is {:?}", B.compress());
 //           println!("B_delta is {:?}", B_delta.compress());
 
            let y_c = a_opening * z_l;
            let P_c = B * a_opening + u_update * y_c;

            assert_eq!(B.compress(), B_delta.compress());

            assert_eq!(P_c.compress(), P.compress());

            return true;  
        } else {
//            println!("======P is {:?}", P.compress());
        }

        let n_prime =  n/2;
        let k_vec_prime = k_vec[1..].to_vec(); //removes the first element

    //1    L_vec[j].compress().append_to_transcript(b"PCS point L", transcript);
    //1   R_vec[j].compress().append_to_transcript(b"PCS point R", transcript);

    //2    let k_vec = transcript.challenge_vector(b"generate challenge k_vec", n.ilog2() as usize);        

       // println!("verifier k_vec {:?}, for gen len {:?}", k_vec, n);
//        println!("verifier B_R_vec_vec[0] {:?}", B_RR_vec_vec[0]);
    
//        println!("^^^^^^j is {:?}", j);
//        println!("^^^^^^k_vec is {:?}", k_vec);
//        println!("^^^^^^B_RR_vec_vec is {:?}", B_RR_vec_vec);
        let mut B_RR_vec_vec_j: Vec<GroupElement> = Vec::new();

        if B_RR_vec_vec.len() > j {
            B_RR_vec_vec_j  = B_RR_vec_vec[j].iter().map(|C| C.decompress().unwrap()).collect();
        }
    //2    let B_delta_updated = B_delta + GroupElement::vartime_multiscalar_mul(&k_vec, &B_R_vec_vec_j);    
    //5    let B_delta_R = B_R + B_R * k_vec[j] + GroupElement::vartime_multiscalar_mul(k_vec[j..].to_vec(), B_RR_vec);    
       // println!("beta_vec is {:?}", beta_vec);
        
    //2    for _index in 0..beta_vec.len() {
    //2        beta_vec[_index] += k_vec[_index];
    //2    }

    //1    B_delta_R_vec[j].compress().append_to_transcript(b"PCS point B_delta_R", transcript);
    
 //   println!("^^^^^^B_R_vec.len() is {:?}", B_R_vec.len());
 //   println!("^^^^^^B_RR_vec_vec_j.len() is {:?}", B_RR_vec_vec_j.len());
 //   println!("^^^^^^k_vec[j+1..].to_vec().len() is {:?}", k_vec[j+1..].to_vec().len());
        let mut B_delta_R = B_R_vec[j] + B_R_vec[j] * k_vec[j]; //+  GroupElement::vartime_multiscalar_mul(k_vec[j..].to_vec(), B_RR_vec_vec_j); 
        if B_RR_vec_vec_j.len() > 0 {
            
            B_delta_R = B_delta_R +  GroupElement::vartime_multiscalar_mul(k_vec[j+1..].to_vec(), B_RR_vec_vec_j); 
        }

        let x = transcript.challenge_scalar(b"generate challenge x");
        //let x = Scalar::one() + Scalar::one() + Scalar::one();
        //let x = Scalar::one() ;
        let x_inv = x.invert().unwrap();

 //       println!("Verfier x is {:?} x_inv is {:?} ", x, x_inv);

       // println!("~~~~~~~~L_vec[j] is {:?} R_vec[j] is {:?} ", L_vec[j].compress(), R_vec[j].compress());

        let P_prime = L_vec[j] * x_inv + P + R_vec[j] * x;


     //   println!("--------------------------------P is {:?}, P_prime is {:?}", P.compress(), P_prime.compress());


        let B_prime = B - B_R_vec[j] + B_R_vec[j] * x_inv;
     //5   let B_prime = B - B_R + B_R * x_inv;
     //2   let B_delta_prime = (B_delta_updated - B_delta_R_vec[j]) + (B_delta_R_vec[j] - B_R_vec_vec_j[0] * beta_vec[0]) * x_inv;
    //5    let B_delta_prime = B_delta - B_delta_R + (B_delta_R - B_R * k_vec[j]) * x_inv;


        let B_delta_prime = B_delta - B_delta_R + (B_delta_R - B_R_vec[j] * k_vec[j]) * x_inv; //wrong, just a place holder

       // println!("z_l is {:?}, z_r is {:?}", z_l, z_r);
        let z_l_prime = z_l + z_r * x_inv;
        let mut z_r_prime = Scalar::zero();

        if (n_prime != 1) {
          //  println!("e_vec.len() is {:?} n.ilog2() is {:?} j is {:?}", e_vec.len(), n.ilog2(), j);
            z_r_prime = z_l_prime * e_vec[(e_vec.len() - 2) as usize-j];
            //println!("e_vec[(n.ilog2() - 1) as usize-j is {:?}", e_vec[(n.ilog2() - 2) as usize-j]);
        }

   //     println!("z_l' is {:?}, z_r' is {:?}", z_l, z_r);

        Self::verifier_verify_recursion(n_prime, u_update, &P_prime, /*2&B_delta_prime,*/ &B_prime, 
             z_l_prime, z_r_prime, y, j+1, &e_vec, &mut beta_vec[1..].to_vec(),
            L_vec, R_vec,/*3  B_delta_R_vec,*/ B_R_vec, B_RR_vec_vec, &B_delta_prime, k_vec, a_opening, transcript) 
    }

}

#[cfg(test)]
mod test {
    //use curve25519_dalek::traits::VartimeMultiscalarMul;

    use crate::groups::transcript;

    use super::*;


    #[test]
    fn test_prove_verify() {
        let mut prover_transcript = Transcript::new(b"unit test");
        let polys = PCSGens::new(256, b"gens_polyeval_verify");
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
    

  //   println!("*****B_R_vec_vec.len() is {:?}", B_RR_vec_vec.len());
        
        let mut verifier_transcript = Transcript::new(b"unit test");

        let (B_delta, B_R, B_RR_vec, k_vec) = PCSEval::verifier_setup(&polys.gens.G);
  
     //   println!("*****B_R is {:?}", B_R);
    //    println!("*****B_R_vec.len() is {:?}", B_R_vec.len());
        let full_B_R_vec = [vec![B_R],  B_R_vec].concat();
        let full_B_RR_vec_vec = [vec![B_RR_vec],  B_RR_vec_vec].concat();
        
        
        PCSEval::verifier_verify(&polys.gens.G, &polys.gens.h, &C, &B, &full_B_R_vec, &full_B_RR_vec_vec, z_1, y, 
            &L_vec, &R_vec, &B_delta, &k_vec, &a_opening, &mut verifier_transcript);

    
    }
/* 
    #[test]
    fn test_prove_2() {
        let mut prover_transcript = Transcript::new(b"unit test");
        let polys = PCSGens::new(2, b"gens_polyeval_verify");
        let mut csprng: OsRng = OsRng;

        let (B, B_R_vec) = PCSEval::prover_setup(&polys.gens.G);

        let x1 = Scalar::one() + Scalar::one() + Scalar::one() + Scalar::one() + Scalar::one();

        let a_1 = Scalar::random(&mut csprng);
        let a_2 = Scalar::random(&mut csprng);

        let g_1 = &polys.gens.G[0]; 
        let g_2 = &polys.gens.G[1]; 

        let u = &polys.gens.h * x1;

        let z_1 = Scalar::one() + Scalar::one();
        let z_2 = z_1 + z_1;

        

        let x = Scalar::one()+ Scalar::one() + Scalar::one();
        let x_inv = x.invert().unwrap();

        let a_vec = vec![a_1, a_2];

        let C = g_1 * a_1 + g_2 * a_2;


        let y = a_1 * z_1 + a_2 * z_2; 

        let (L_vec, R_vec, B_delta_R_vec, B_R_vec_vec, a_opening, z_opening) 
            = PCSEval::prover_prove(&polys.gens.G, &polys.gens.h, &C, &B, &B_R_vec, z_1, y, &a_vec, &mut prover_transcript);
    
        let mut verifier_transcript = Transcript::new(b"unit test");
        PCSEval::verifier_verify(&polys.gens.G, &polys.gens.h, &C, &B, z_1, y, 
            &L_vec, &R_vec, &B_delta_R_vec, &B_R_vec_vec, &a_opening, &mut verifier_transcript);

        let L = g_2 * a_1 + u * (a_1 * z_2);
        let R = g_1 * a_2 + u * (a_2 * z_1);

        let P = C + u * y;
        let P_prime = L * x_inv + P + R * x;

        let B = g_1 + g_2 * x_inv;
        let a_prime = a_1 + a_2 * x;
/* 
        println!("------- L = {:?}, R = {:?}", L.compress(), R.compress());
        
        println!("------- P = {:?} P_prime - {:?}", P.compress(), P_prime.compress());
        
        println!("------- B = {:?} a_prime - {:?}", B.compress(), a_prime);
   */     
    
    }
*/


    #[test]
    fn test_compute_L_R_2() {

        let mut csprng: OsRng = OsRng;
        let x = Scalar::random(&mut csprng);
        
        let a_1 = Scalar::random(&mut csprng);
        let a_2 = Scalar::random(&mut csprng);

        println!("a_1 is {:?}", a_1);
        println!("a_2 is {:?}", a_2);

        let z_1 = Scalar::random(&mut csprng);
        let z_2 = Scalar::random(&mut csprng);


        let a_vec = vec![a_1, a_2];
        let z_vec = vec![z_1, z_2];

        let polys = PCSGens::new(2, b"gens_polyeval_verify");

        let gens = polys.gens.G;
        let u = polys.gens.h;
        let n_prime =  gens.len()/2;
        let gens_L = vec![gens[1]];
        let gens_R = vec![gens[0]];

        let (L_c,R_c) = PCSEval::compute_L_R(&gens_L, &gens_R, &polys.gens.h, &a_vec, &z_vec);

        let c_L = a_1 * z_2 ;
        let L =  gens[1] * a_1 + polys.gens.h * c_L;

        let c_R = a_2 * z_1 ;
        let R =  gens[0] * a_2 + polys.gens.h * c_R;

        assert_eq!(L_c.compress(), L.compress());
        assert_eq!(R_c.compress(), R.compress());

        let y = a_1 * z_1 + a_2 * z_2;
        let P = gens[0] * a_1 + gens[1] * a_2 + u * y;
        
        let x = Scalar::one() + Scalar::one() + Scalar::one();
        let x_inv = x.invert().unwrap();

        let P_prime = L * x_inv + P + R * x;

        let g_prime = gens[0] + gens[1] * x_inv;
        let a_prime = a_1 + a_2 * x;
        let z_prime = z_1 + z_2 * x_inv;

        let P_c = g_prime * a_prime + u * (a_prime * z_prime);

        println!("P_c is {:?}", P_c.compress());

        println!("P_prime is {:?}", P_prime.compress());




    }

    #[test]
    fn test_compute_L_R() {

        let mut csprng: OsRng = OsRng;
        let x = Scalar::random(&mut csprng);
        
        let a_1 = Scalar::random(&mut csprng);
        let a_2 = Scalar::random(&mut csprng);
        let a_3 = Scalar::random(&mut csprng);
        let a_4 = Scalar::random(&mut csprng);

        let z_1 = Scalar::random(&mut csprng);
        let z_2 = Scalar::random(&mut csprng);
        let z_3 = Scalar::random(&mut csprng);
        let z_4 = Scalar::random(&mut csprng);

        let a_vec = vec![a_1, a_2, a_3, a_4];
        let z_vec = vec![z_1, z_2, z_3, z_4];

        let polys = PCSGens::new(4, b"gens_polyeval_verify");

        let gens = polys.gens.G;
        let n_prime =  gens.len()/2;
        let gens_L = gens[..n_prime].to_vec();
        let gens_R = gens[n_prime..].to_vec();

        let (L_c,R_c) = PCSEval::compute_L_R(&gens_L, &gens_R, &polys.gens.h, &a_vec, &z_vec);

        let c_L = a_1 * z_3 + a_2 * z_4;
        let L =  gens[2] * a_1 + gens[3] * a_2 + polys.gens.h * c_L;

        let c_R = a_3 * z_1 + a_4 * z_2;
        let R =  gens[0] * a_3 + gens[1] * a_4 + polys.gens.h * c_R;

        assert_eq!(L_c.compress(), L.compress());
        assert_eq!(R_c.compress(), R.compress());

    }


    #[test]
    fn test_compute_primes_vec() {

        let mut csprng: OsRng = OsRng;
        let x = Scalar::random(&mut csprng);
        let one_1 = Scalar::random(&mut csprng);
        let one_2 = Scalar::random(&mut csprng);
        let one_3 = Scalar::random(&mut csprng);
        let one_4 = Scalar::random(&mut csprng);
        let one_5 = Scalar::random(&mut csprng);
        let one_6 = Scalar::random(&mut csprng);
        let one_7 = Scalar::random(&mut csprng);
        let one_8 = Scalar::random(&mut csprng);


        let a_vec = vec![one_1, one_2, one_3, one_4, one_5, one_6, one_7, one_8];

        let a_prime_vec = PCSEval::compute_prime_vec(&a_vec, &x);

        let a_prime_1 = a_vec[0] + a_vec[4] * x;
        let a_prime_2 = a_vec[1] + a_vec[5] * x;
        let a_prime_3 = a_vec[2] + a_vec[6] * x;
        let a_prime_4 = a_vec[3] + a_vec[7] * x;

        assert_eq!(a_prime_vec[0], a_prime_1);
        assert_eq!(a_prime_vec[1], a_prime_2);
        assert_eq!(a_prime_vec[2], a_prime_3);
        assert_eq!(a_prime_vec[3], a_prime_4);
    }

    #[test]
    fn test_compute_g_primes_vec() {

        let polys = PCSGens::new(8, b"gens_polyeval_verify");

        let mut csprng: OsRng = OsRng;
        let x = Scalar::random(&mut csprng);

        let g_prime_vec = PCSEval::compute_g_prime_vec(&polys.gens.G, &x);

        let g_prime_1 = polys.gens.G[0] + polys.gens.G[4] * x;
        let g_prime_2 = polys.gens.G[1] + polys.gens.G[5] * x;
        let g_prime_3 = polys.gens.G[2] + polys.gens.G[6] * x;
        let g_prime_4 = polys.gens.G[3] + polys.gens.G[7] * x;

        assert_eq!(g_prime_vec[0].compress(), g_prime_1.compress());
        assert_eq!(g_prime_vec[1].compress(), g_prime_2.compress());
        assert_eq!(g_prime_vec[2].compress(), g_prime_3.compress());
        assert_eq!(g_prime_vec[3].compress(), g_prime_4.compress());
    }

    #[test]
    fn test_compute_halves() {

        let polys = PCSGens::new(8, b"gens_polyeval_verify");
        
        let B_Rs = PCSEval::compute_halves(&polys.gens.G);
        let B_R_one = polys.gens.G[4] + polys.gens.G[5] + polys.gens.G[6] +polys.gens.G[7];
        let B_R_two = polys.gens.G[2] + polys.gens.G[3] + polys.gens.G[6] +polys.gens.G[7];
        let B_R_three = polys.gens.G[1] + polys.gens.G[3] + polys.gens.G[5] +polys.gens.G[7];

        assert_eq!(B_R_one.compress(), B_Rs[0].compress());
        assert_eq!(B_R_two.compress(), B_Rs[1].compress());
        assert_eq!(B_R_three.compress(), B_Rs[2].compress());
    }

    #[test]
    fn test_computing_z() {

        let x = Scalar::from(2u64);
        let e_vec = PCSEval::compute_e_vec(x, 32);
        println!("e_vec is {:?}", e_vec);

        let one = x * x.invert().unwrap();

        println!("one is {:?} from {:?}", one, x.invert().unwrap());

    }
    #[test]
    fn test_computing_exponents() {

        let two = Scalar::one() + Scalar::one();
        let three = two + Scalar::one();
        let four = two + two;

        let alpha_vec = vec![Scalar::one();16];
        let k_vec = vec![Scalar::one(), two, three, four];

        let alpha_prime_vec = PCSEval::compute_exponents(&alpha_vec, &k_vec);
        

        println!("alpha_prime_vec is {:?} ", alpha_prime_vec);

    }

}