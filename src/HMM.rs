// lib.rs

use crate::VDJmodeler::SequenceModel;
use crate::VDJmodeler::HMMmodel;
 use crate::VDJmodeler::HMMcollector;

pub struct HMMState {
    pub match_emission: Vec<Vec<f64>>,    // Probabilities for A, G, C, T, and `.`
}

impl HMMState {

    /*
    // Method to get the index of the maximum emission probability
    pub fn max_emission(&self) -> f64 {
        let mut max_value = self.match_emission[0];

        // Find the maximum in match_emission
        for &value in &self.match_emission {
            if value > max_value {
                max_value = value;
            }
        }

        max_value
    }
    // Method to get the index of the maximum emission probability
    pub fn max_emission_index(&self) -> usize {
        let mut max_index = 0;
        let mut max_value = self.match_emission[0];

        // Find the maximum in match_emission
        for (i, &value) in self.match_emission.iter().enumerate() {
            if value > max_value {
                max_value = value;
                max_index = i;
            }
        }

        // Check if there is a higher value in insertion_emission
        for (i, &value) in self.insertion_emission.iter().enumerate() {
            if value > max_value {
                max_value = value;
                max_index = i; // Offset by 5 to differentiate from match_emission indices
            }
        }

        max_index
    }

    // Method to get the index of the maximum emission probability
    pub fn emission_index_for(&self, value: &f64) -> usize {


        // Find the maximum in match_emission
        for (i, &mine) in self.match_emission.iter().enumerate() {
            if mine == *value {
                return i
            }
        }
        panic!("The value {value} is not in my list!");
    }
    */
    fn adjust_zero_probabilities(probabilities: &mut [f64], min_prob: f64) {
        for prob in probabilities.iter_mut() {
            if *prob == 0.0 {
                *prob = min_prob;
            }
        }
    }
    fn normalize_probabilities(probabilities: &mut [f64]) {
        let total: f64 = probabilities.iter().sum();
        if total > 0.0 {
            for prob in probabilities.iter_mut() {
                *prob /= total;
            }
        }
    }

    // Normalize match emission probabilities based on a list of HMMcollectors
    pub fn from_collectors( collectors: &[HMMcollector]) -> Self {
        //pub match_emission: Vec<[f64; 5]>,    // Probabilities for A, G, C, T, and `.`
        //pub insertion_emission: Vec<[f64; 5]>,
        let mut total_counts = vec![0; 5];
        
        // Aggregate counts from all collectors
        for collector in collectors {
            for (i, &count) in collector.states.iter().enumerate() {
                total_counts[i] += count;
            }
        }
        let total = total_counts.iter().sum::<usize>() as f64;

        let mut match_emission = vec![vec![0.0; 5]; collectors.len() ];
        let min_prob = 0.0001;
        for (i, collector) in collectors.iter().enumerate() {
            for (j, &count) in collector.states.iter().enumerate() {
                match_emission[i][j] = count as f64 /total;
            }
            Self::adjust_zero_probabilities( &mut match_emission[i], min_prob );
            Self::normalize_probabilities( &mut match_emission[i]);
        }

        Self{
            match_emission,
        }
    }

}

pub struct HMM {
    /// the different HMM emission states for the different names at each position
    pub states: Vec<HMMState>,
    /// the fixed transition matrix - never switch from IGH to e.g. TRA
    pub transition_matrix: Vec<Vec<f64>>,
    /// the different sequence names
    pub names: Vec<SequenceModel>,
}

impl HMM {
    // Create an HMM from the given sequence models
    pub fn from_sequence_models(models: Vec<HMMmodel>) -> Self {
        // Assume all models have the same length for simplicity
        let sequence_length = models.iter().map(|m| m.collector.len()).max().unwrap_or(0);

        // Define a fixed transition matrix
        // Example: Suppose we have `num_states` states (one for each model)
        let num_states = models.len(); 
        let low_prob = 0.001;
        let diagonal_prob = 1.0 - (low_prob * (num_states as f64 - 1.0));

        // Create transition matrix with the calculated probabilities
        let mut transition_matrix = vec![vec![low_prob; num_states]; num_states];
        for i in 0..num_states {
            transition_matrix[i][i] = diagonal_prob;
        }

        // Create HMM states from models
        let states: Vec<HMMState> = (0..sequence_length)
            .map(|i| {
                let collectors: Vec<HMMcollector> = models.iter()
                    .filter_map(|model| model.collector.get(i).cloned())
                    .collect();
                HMMState::from_collectors(&collectors)
            }).collect();

        // Create HMM instance
        HMM { 
            states,
            transition_matrix,
            names: models.iter().map(|m| m.name.clone()).collect(), // Clone the names
        }
    }

    /// translate the sequence into the correct position in the data
    pub fn char2pos( seq: char) -> usize{
        match seq {
            'A' => 0,
            'G' => 1,
            'C' => 2,
            'T' => 3,
            '.' => 4,
            _ => unreachable!()
        }
    }

    pub fn states( &self ) -> &Vec<HMMState>{
        &self.states
    }

    /// Forward algorithm
    /// The only thing I really need from this as I 'only' want to check if any of the sequences
    /// would be of a VDJ recombination evet.
    pub fn forward_algorithm(&self, sequence: &str) -> Vec<(String, f64)> {
        let num_states = self.states.len();
        let sequence_length = sequence.len();

        let mut alpha = vec![vec![0.0; num_states]; sequence_length];

        // Initialize the alpha values for the first position
        for state in 0..num_states {
            let emission_prob = &self.states[0].match_emission[Self::char2pos(sequence.chars().nth(0).unwrap()) ];
            alpha[0][state] = emission_prob[state];
        }

        // Recursively compute the alpha values for the rest of the sequence
        for t in 1..sequence_length {
            for j in 0..num_states {
                alpha[t][j] = 0.0;
                for i in 0..num_states {
                    alpha[t][j] += alpha[t - 1][i] * self.transition_matrix[i][j];
                }
                let emission_prob = &self.states[t].match_emission[Self::char2pos(sequence.chars().nth(t).unwrap()) ];
                alpha[t][j] *= emission_prob[j];
            }
        }

        // Get the final probabilities from the last position in the alpha matrix
        let final_probabilities: Vec<f64> = alpha[sequence_length - 1].clone();

        // Compute the final probability
        let total_prob: f64 = final_probabilities.iter().sum();

        // Pair the names with their corresponding probabilities
        let mut result: Vec<(String, f64)> = self.names.iter()
            .zip(final_probabilities.iter())
            .map(|(name, &prob)| (name.name(), prob))
            .collect();
        // Add the total probability as the last entry
        result.push(("Total".to_string(), total_prob));

        result
    }

    /*
    // Backward algorithm
    pub fn backward_algorithm(&self, sequence: &[usize]) -> f64 {
        let num_states = self.states.len();
        let sequence_length = sequence.len();

        let mut beta = vec![vec![0.0; num_states]; sequence_length];

        // Initialize the beta values for the last position
        for state in 0..num_states {
            beta[sequence_length - 1][state] = 1.0;
        }

        // Recursively compute the beta values for the rest of the sequence
        for t in (0..sequence_length - 1).rev() {
            for i in 0..num_states {
                beta[t][i] = 0.0;
                for j in 0..num_states {
                    let emission_prob = self.states[t + 1].match_emission[sequence[t + 1]];
                    beta[t][i] += self.transition_matrix[i][j] * emission_prob[j] * beta[t + 1][j];
                }
            }
        }

        // Compute the final probability
        let mut prob = 0.0;
        for state in 0..num_states {
            let emission_prob = self.states[0].match_emission[sequence[0]];
            prob += emission_prob[state] * beta[0][state];
        }
        prob
    }

    // Final function to get the most likely sequence of hidden states
    pub fn most_likely_states(&self, sequence: &[usize]) -> Vec<usize> {
        let num_states = self.states.len();
        let sequence_length = sequence.len();

        let mut delta = vec![vec![0.0; num_states]; sequence_length];
        let mut psi = vec![vec![0; num_states]; sequence_length];

        // Initialize delta and psi for the first position
        for state in 0..num_states {
            let emission_prob = self.states[0].match_emission[sequence[0]];
            delta[0][state] = emission_prob[state];
            psi[0][state] = 0;
        }

        // Recursively compute delta and psi for the rest of the sequence
        for t in 1..sequence_length {
            for j in 0..num_states {
                let mut max_delta = 0.0;
                let mut best_state = 0;
                for i in 0..num_states {
                    let prob = delta[t - 1][i] * self.transition_matrix[i][j];
                    if prob > max_delta {
                        max_delta = prob;
                        best_state = i;
                    }
                }
                delta[t][j] = max_delta * self.states[t].match_emission[sequence[t]][j];
                psi[t][j] = best_state;
            }
        }

        // Backtrack to find the most likely sequence of states
        let mut most_likely_states = vec![0; sequence_length];
        let mut last_state = delta[sequence_length - 1].iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap()).unwrap().0;
        most_likely_states[sequence_length - 1] = last_state;

        for t in (0..sequence_length - 1).rev() {
            last_state = psi[t + 1][last_state];
            most_likely_states[t] = last_state;
        }

        most_likely_states
    }
    */
    // Placeholder for functions like forward_algorithm, viterbi, etc.
}

