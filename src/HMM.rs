// lib.rs

use crate::VDJmodeler::SequenceModel;
use crate::VDJmodeler::HMMmodel;
use crate::VDJmodeler::HMMcollector;
use std::collections::HashSet;
use std::collections::HashMap;

use std::f64;

pub struct HMMState {
    pub match_emission: Vec<Vec<f64>>,    // Probabilities for A, G, C, T, and `.`
}

impl HMMState {

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

    pub fn prob_for_pos(&self, pos:usize) -> Vec<f64> {
        let mut ret = Vec::<f64>::with_capacity(self.match_emission.len());
        for emissions in &self.match_emission {
            ret.push (*emissions.get(pos).unwrap_or( &0.00001 ));
        }
        ret
    }

    pub fn len(&self) -> usize {
        self.match_emission.len()
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
    pub fn char2pos( seq: u8) -> Option<usize>{
        match seq.to_ascii_uppercase() {
            b'A' => Some(0),
            b'G' => Some(1),
            b'C' => Some(2),
            b'T' => Some(3),
            other => {
                //eprintln!("Sorry I can not decode this char {}",other as char );
                None
            },
        }
    }

    pub fn iupac_char2pos(seq: u8) -> Vec<usize> {
    match seq.to_ascii_uppercase() {
        b'A' => vec![0],        // A or a
        b'G' => vec![1],        // G or g
        b'C' => vec![2],        // C or c
        b'T' => vec![3],        // T or t
        b'.' => vec![0, 1, 2, 3],        // . (gap)
        b'N' => vec![0, 1, 2, 3], // N or n
        b'R' => vec![0, 1],     // R or r
        b'Y' => vec![2, 3],     // Y or y
        b'W' => vec![0, 3],     // W or w
        b'S' => vec![1, 2],     // S or s
        b'K' => vec![1, 3],     // K or k
        b'M' => vec![0, 2],     // M or m
        b'B' => vec![1, 2, 3],  // B or b
        b'D' => vec![0, 1, 3],  // D or d
        b'H' => vec![0, 2, 3],  // H or h
        b'V' => vec![0, 1, 2],  // V or v
        other => panic!("Sorry, I cannot decode this char: {}", other as char),
    }
}

    pub fn states( &self ) -> &Vec<HMMState>{
        &self.states
    }

    fn try_start_at(&self, sequence: &[u8], pos: usize) -> Option<Vec<f64>> {
        // Initialize the return vector with zeros
        let mut ret = vec![0.0; self.states[0].match_emission.len()];

        // Iterate over the sequence
        for t in 0..sequence.len() {
            let seq_id = match HMM::char2pos(sequence[t]) {
                Some(id) => id,
                None => {
                    // Return a vector of zeros if there's an invalid character
                    return None
                }
            };

            // Accumulate probabilities into ret
            let prob_for_pos = self.states[t+pos].prob_for_pos(seq_id);

            // Ensure the lengths match before adding
            if prob_for_pos.len() != ret.len() {
                //eprintln!("Warning: Length mismatch between probability vector and return vector at position {} with start {pos}", t+pos );
                return None
            }

            // Perform element-wise addition
            for (ret_val, prob_val) in ret.iter_mut().zip(prob_for_pos.iter()) {
                *ret_val += prob_val;
            }
        }

        Some(ret)
    }

    pub fn find_probable_start(&self, sequence: &[u8]) -> Vec<(usize, f64)> {
        let num_states = self.states[0].len();
        let mut stats_vec = vec![0.0; num_states];
        let mut pos_vec = vec![0; num_states];

        // Ensure that sequence length is not greater than the number of positions available
        let max_start_pos = self.states.len().saturating_sub(sequence.len());

        for t in 0..max_start_pos {
            let current_probs = match self.try_start_at(sequence, t){
                Some(ret) => ret,
                None => {
                    //out of data in the model;
                    break;
                }
            };

            for (id, value) in current_probs.iter().enumerate() {
                if stats_vec[id] < *value {
                    stats_vec[id] = *value;
                    pos_vec[id] = t;
                }
            }
        }

        // Combine pos_vec and stats_vec into a Vec<(usize, f64)>

        for i in 0..stats_vec.len(){
            stats_vec[i] =  stats_vec[i] / sequence.len() as f64;
        }
        pos_vec.into_iter()
            .zip(stats_vec.into_iter())
            .collect()
        
    }

    /// Forward algorithm
    /// The only thing I really need from this as I 'only' want to check if any of the sequences
    /// would be of a VDJ recombination evet.
    pub fn forward_algorithm(&self, sequence: &[u8]) -> Option< Vec<(String, f64)> > {
        

        let probable_start_values = self.find_probable_start( sequence );

        let mut start_values:HashSet<usize> = HashSet::new();
        for  (pos, stat ) in &probable_start_values{
            if stat > &0.3 {
                start_values.insert( *pos);
            }
        }

        //let start_values:Vec<usize> = probable_start_values.iter().filter_map(|x| { if x.1 > 0.3{ Some( x.0 ) }else { None} }).collect();

        if start_values.len() == 0 {
            return None
        }
        //println!("Using these probable_start_values {probable_start_values:?}\nI find these probable start positions: {start_values:?}");

        let mut data: Vec<Option<Vec<(String, f64)>>> = vec![];
        for start in start_values {
            let this = self.forward_algorithm_pos( sequence, start );
            data.push( this );
        }

        Some(HMM::collapse_to_max( data ))
    }

    fn collapse_to_max(data: Vec<Option<Vec<(String, f64)>>> ) -> Vec<(String, f64)> {
        let mut max_values: HashMap<String, f64> = HashMap::new();

        for item in data {
            if let Some(stats) = item {
                for (key, value) in stats {
                    // Update the max value for each key
                    max_values.entry(key)
                        .and_modify(|e| *e = (*e).max(value))
                        .or_insert(value);
                }
            }
        }

        // Convert the HashMap back to a Vec<(String, f64)>
        max_values.into_iter()
            .map(|(key, value)| (key.to_string(), value))
            .collect()
    }

    pub fn forward_algorithm_pos(&self, sequence: &[u8], start:usize) -> Option<Vec<(String, f64)>> {

        let num_states = self.states[0].len();
        let mut sequence_length = sequence.len();

        // do we have enough info in the model to do this:
        let this_end = self.states.len().min(start + sequence_length );

        if start > this_end {
            return None
        }

        let mut alpha = vec![vec![-f64::INFINITY; num_states]; sequence_length];
        // Initialize the alpha values for the first position

        let first_pos = match Self::char2pos(sequence[0]){
            Some(pos) => pos,
            None => {
                //eprintln!("Bad sequence:")
                return None
            }
        };

        for (state, emission_prob) in self.states[start].prob_for_pos( first_pos ).iter().enumerate(){
             alpha[0][state] = emission_prob.ln();
        }


        // Recursively compute the alpha values for the rest of the sequence
        'main: for t in start+1..this_end {
            let current_pos = match Self::char2pos(sequence[t-start]){
                Some(pos) => pos,
                None => return None,
            };
            for (state, emission_prob) in self.states[t].prob_for_pos( current_pos ).iter().enumerate(){
                let mut max_log = -f64::INFINITY;
                for i in 0..num_states {
                    max_log = max_log.max(alpha[t - 1-start][i] + self.transition_matrix[i][state].ln());
                }
                alpha[t-start][state] = max_log + emission_prob.ln();
            }

        }
        sequence_length = this_end - start;
        // Compute the final probability in log-space
        let final_probabilities: Vec<f64> = alpha[sequence_length - 1].clone();
        let max_final_prob = final_probabilities.iter().cloned().fold(-f64::INFINITY, f64::max);
        let total_prob = (max_final_prob).exp() ;

        // Pair the names with their corresponding probabilities
        let mut result: Vec<(String, f64)> = self.names.iter()
            .zip(final_probabilities.iter())
            .map(|(name, &prob)| (name.name(), (prob - max_final_prob).exp()))
            .collect();
        // Add the total probability as the last entry
        //result.push(("Total".to_string(), total_prob));

        Some(result)
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

