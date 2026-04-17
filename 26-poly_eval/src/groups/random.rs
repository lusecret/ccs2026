use super::super::scalar::Scalar;
use super::transcript::ProofTranscript;
use merlin::Transcript;
use rand::rngs::OsRng;

pub struct RandomTape {
  tape: Transcript,
}

impl RandomTape {
  pub fn new(name: &'static [u8]) -> Self {
    let tape = {
      let mut csprng: OsRng = OsRng;
      let mut tape = Transcript::new(name);
      tape.append_scalar(b"init_randomness", &Scalar::random(&mut csprng));
      tape
    };
    Self { tape }
  }

  pub fn new_ME(name: &'static [u8]) -> Self {
    let tape = {
      let mut tape = Transcript::new(name);
      tape.append_scalar(b"init_randomness", &Scalar::from_bytes(b"------initial randomness from OS").unwrap());
      tape
    };
    Self { tape }
  }

  pub fn random_scalar(&mut self, label: &'static [u8]) -> Scalar {
    self.tape.challenge_scalar(label)
  }

  pub fn random_scalar_by_size(&mut self, label: &'static [u8], size: usize) -> Scalar {
    self.tape.challenge_by_size(label, size)
  }

  pub fn random_vector(&mut self, label: &'static [u8], len: usize) -> Vec<Scalar> {
    self.tape.challenge_vector(label, len)
  }

  pub fn random_vector_q(&mut self, label: &'static [u8], len: usize) -> Vec<Scalar> {
    self.tape.challenge_vector_mont(label, len)
  }

  pub fn random_vector_by_size(&mut self, label: &'static [u8], len: usize, size: usize) -> Vec<Scalar> {
    self.tape.challenge_vector_by_size(label, len, size)
  }
}
