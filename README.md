# ccs2026
code used to benchmark protocol for CCS2026 submission titled: "Efficient Transparent Polynomial Commitment Scheme Evaluated with Unknown Constraint"

##Quick Start

1.	Download rust
   
3.	Use cargo bench command to run benchmark testing
   
##Benchmarking Options

1.	Go to file: 26-poly_eval/benches/poly_eval_bench.rs
   
3.	Uncomment the section if you want to benchmark (prove or verify)
targets = poly_eval_prove_benchmark
//targets = poly_eval_verify_benchmark

4.	In the default setting, we set the polynomial size to 2^20, you can reset it to other test size by resetting the value inside the brackets [ ] 
for &s in [20].iter() { //valid values are 10, 12, 14, 16, 18, 20

<img width="468" height="305" alt="image" src="https://github.com/user-attachments/assets/742f93a4-6042-4c9d-b7f0-cb5ecf63b4ef" />
