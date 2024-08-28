// lib.rs


pub struct HMMState {
    pub match_emission: [f64; 5],    // Probabilities for A, G, C, T, and `.`
    pub insertion_emission: [f64; 5],
    pub transition_probabilities: [f64; 3], // Example: Match->Match, Match->Insertion, etc.
}

impl HMMState {
    fn new() -> Self {
        HMMState {
            match_emission: [0.2, 0.2, 0.2, 0.2, 0.2], // Uniform probabilities
            insertion_emission: [0.25, 0.25, 0.25, 0.25, 0.0], // Example probabilities
            transition_probabilities: [0.9, 0.05, 0.05], // Example transition probabilities
        }
    }

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
}

pub struct HMM {
    pub states: Vec<HMMState>,
}

impl HMM {
    pub fn new(num_states: usize) -> Self {
        let states = (0..num_states).map(|_| HMMState::new()).collect();
        HMM { states }
    }

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

    pub fn forward_algorithm(&self, sequence: &str) -> Vec<Vec<f64>> {
        let num_states = self.states.len();
        let seq_len = sequence.len();
        let mut alpha = vec![vec![0.0; num_states]; seq_len];

        // Initialize
        let first_char = sequence.chars().next().unwrap_or('.');
        let first_char_index = HMM::char2pos( first_char);
        for (i, state) in self.states.iter().enumerate() {
            alpha[0][i] = state.match_emission[first_char_index] * 1.0; // Assuming initial probability is 1.0
        }

        // Recurrence
        for t in 1..seq_len {
            let char_index = HMM::char2pos(  sequence.chars().nth(t).unwrap_or('.') );
            for j in 0..num_states {
                alpha[t][j] = self.states.iter().enumerate().map(|(i, state)| {
                    alpha[t - 1][i] * state.transition_probabilities[0] // Simplified transition
                }).sum::<f64>() * self.states[j].match_emission[char_index];
            }
        }

        alpha
    }

    pub fn backward_algorithm(&self, sequence: &str) -> Vec<Vec<f64>> {
        let num_states = self.states.len();
        let seq_len = sequence.len();
        let mut beta = vec![vec![0.0; num_states]; seq_len];

        // Initialize
        for i in 0..num_states {
            beta[seq_len - 1][i] = 1.0; // End state probability is 1.0
        }

        // Recurrence
        for t in (0..(seq_len - 1)).rev() {
            let char_index = HMM::char2pos(  sequence.chars().nth(t + 1).unwrap_or('.') );
            for i in 0..num_states {
                beta[t][i] = self.states.iter().enumerate().map(|(j, state)| {
                    beta[t + 1][j] * state.transition_probabilities[0] * self.states[j].match_emission[char_index]
                }).sum::<f64>();
            }
        }

        beta
    }

    pub fn baum_welch(&mut self, sequences: Vec<&str>, num_iterations: usize) {
        let num_states = self.states.len();
        
        for _ in 0..num_iterations {
            let mut a_count = vec![vec![0.0; num_states]; num_states];
            let mut b_count = vec![vec![0.0; 5]; num_states]; // Counts for emissions
            
            for sequence in &sequences {
                let alpha = self.forward_algorithm(sequence);
                let beta = self.backward_algorithm(sequence);

                let seq_len = sequence.len();
                
                // E-step: Calculate expected counts
                for t in 0..(seq_len - 1) {
                    let char_index = HMM::char2pos(  sequence.chars().nth(t).unwrap_or('.') );

                    for i in 0..num_states {
                        for j in 0..num_states {
                            a_count[i][j] += alpha[t][i] * self.states[i].transition_probabilities[0] * self.states[j].match_emission[char_index] * beta[t + 1][j];
                        }
                        b_count[i][char_index] += alpha[t][i] * beta[t][i];
                    }
                }
            }

            // M-step: Update HMM parameters
            for i in 0..num_states {
                let denom = a_count[i].iter().sum::<f64>();
                for j in 0..num_states {
                    self.states[i].transition_probabilities[j] = a_count[i][j] / denom;
                }

                let b_denom = b_count[i].iter().sum::<f64>();
                for k in 0..5 {
                    self.states[i].match_emission[k] = b_count[i][k] / b_denom;
                }
            }
        }
    }

    pub fn viterbi(&self, sequence: &str) -> (Vec<usize>, f64) {
        let num_states = self.states.len();
        let seq_len = sequence.len();
        let mut viterbi = vec![vec![0.0; num_states]; seq_len];
        let mut backpointer = vec![vec![0; num_states]; seq_len];

        // Initialization
        let first_char_index = HMM::char2pos(sequence.chars().next().unwrap_or('.'));
        for (i, state) in self.states.iter().enumerate() {
            viterbi[0][i] = state.match_emission[first_char_index] * 1.0; // Assume initial probability is 1.0
            backpointer[0][i] = 0;
        }

        // Recurrence
        for t in 1..seq_len {
            let char_index = HMM::char2pos(sequence.chars().nth(t).unwrap_or('.'));
            for j in 0..num_states {
                let (max_prob, best_state) = (0..num_states).map(|i| {
                    let prob = viterbi[t - 1][i] * self.states[i].transition_probabilities[j]; // Use correct transition probability
                    (prob, i)
                }).max_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)).unwrap();
                
                viterbi[t][j] = max_prob * self.states[j].match_emission[char_index];
                backpointer[t][j] = best_state;
            }
        }

        // Traceback
        let (last_state, max_prob) = viterbi[seq_len - 1].iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal)).unwrap();
        let mut best_path = vec![0; seq_len];
        let mut state = last_state;

        for t in (0..seq_len).rev() {
            best_path[t] = state; // Store state index
            state = backpointer[t][state];
        }

        (best_path, *max_prob as f64)
    }

    // Placeholder for functions like forward_algorithm, viterbi, etc.
}

