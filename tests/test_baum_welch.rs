// test_baum_welch.rs

use hmm_mapper::HMM::{HMM, HMMState};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_baum_welch() {
        // Define a simple HMM with 2 states
        let mut hmm = HMM {
            states: vec![
                HMMState {
                    match_emission: [0.6, 0.2, 0.1, 0.05, 0.05], // A, G, C, T, .
                    insertion_emission: [0.25; 5],
                    transition_probabilities: [0.7, 0.3, 0.0], // From state 0 to 0, 1, or insert
                },
                HMMState {
                    match_emission: [0.1, 0.3, 0.4, 0.1, 0.1],
                    insertion_emission: [0.25; 5],
                    transition_probabilities: [0.4, 0.6, 0.0], // From state 1 to 0, 1, or insert
                }
            ]
        };

        // Sample sequences (each character corresponds to A, G, C, T, .)
        let sequences = vec!["AGCT", "GCTA"];

        // Train the HMM using the Baum-Welch algorithm
        hmm.baum_welch(sequences, 10);

        // Assert updated transition probabilities
        let state0 = &hmm.states()[0];
        let state1 = &hmm.states()[1];

        // Example assertions (these should be adapted based on expected results)
        assert!(state0.transition_probabilities.iter().all(|&p| p >= 0.0 && p <= 1.0), "Transition probabilities should be between 0 and 1");
        assert!(state1.transition_probabilities.iter().all(|&p| p >= 0.0 && p <= 1.0), "Transition probabilities should be between 0 and 1");

        // Emission probabilities should also be between 0 and 1
        assert!(state0.match_emission.iter().all(|&p| p >= 0.0 && p <= 1.0), "Emission probabilities should be between 0 and 1");
        assert!(state1.match_emission.iter().all(|&p| p >= 0.0 && p <= 1.0), "Emission probabilities should be between 0 and 1");

        // Additional assertions based on expected parameter updates could be added here
    }
}