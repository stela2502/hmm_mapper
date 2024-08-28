
// main.rs
use hmm_mapper::HMM::{HMM}; // Import from lib.rs


fn main() {
    let hmm = HMM::new(10); // Initialize HMM with 10 states as an example
    println!("Initialized HMM with {} states.", hmm.states().len());
}

