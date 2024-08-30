//VDHmodeler.rs

use crate::HMM::HMM;
//use crate::fasta_reader::{FastaRecord, FastaReader};
use needletail::parse_fastx_file;
use std::collections::HashSet;

pub enum Chain{
	V,
	D,
	J,
}

impl Chain{
	pub fn id(&self)->usize {
		match self{
			Chain::V => 0,
			Chain::D => 1,
			Chain::J => 2,
		}
	}

	pub fn starts_at(&self, mode:&str, data:&Vec<usize>) -> usize{
		match mode{
			"HeavyChain" => {
				match self{
					Chain::V => 0,
					Chain::D => data[0],
					Chain::J => data[1],
				}
			},
			"LightChain" => {
				match self{
					Chain::V => 0,
					Chain::D => panic!("A light chain has no D segement!"),
					Chain::J => data[0],
					_ => unreachable!(),
				}
			},
			_ => unreachable!(),
		}
	}
}


#[derive(Eq, Hash, PartialEq, Debug, Clone)]
pub enum SequenceModel{
	IGH,
	IGL,
	IGK,
	TRA,
	TRB,
	TRG,
	TRD,
}

impl SequenceModel {
	pub fn id(&self) -> usize{
		match self {
			SequenceModel::IGH => 0,
			SequenceModel::IGL => 1,
			SequenceModel::IGK => 2,
			SequenceModel::TRA => 3, 
			SequenceModel::TRB => 4,
			SequenceModel::TRG => 5,
			SequenceModel::TRD => 6,
			_ => unreachable!(),
		}
	}

	pub fn length() -> usize{
		7
	}

	pub fn from_index(index: usize) -> Option<Self> {
        match index {
            0 => Some(SequenceModel::IGH), // Example, assuming HeavyChain::V
            1 => Some(SequenceModel::IGL), // Example, assuming LightChain::V
            2 => Some(SequenceModel::IGK), // Example, assuming LightChain::V
            3 => Some(SequenceModel::TRA), // Example, assuming HeavyChain::V
            4 => Some(SequenceModel::TRB), // Example, assuming LightChain::V
            5 => Some(SequenceModel::TRG), // Example, assuming HeavyChain::V
            6 => Some(SequenceModel::TRD), // Example, assuming LightChain::V
            _ => None,
        }
    }


	pub fn has_data(&self, data:&Vec<usize> ) -> bool{
		match self {
            SequenceModel::IGH| SequenceModel::TRA | SequenceModel::TRG => {
                // Check for HeavyChain types
                data[0] != 0 && data[1] != 0 && data[2] != 0
            },
            SequenceModel::IGL | SequenceModel::IGK | SequenceModel::TRB | SequenceModel::TRD => {
                // Check for LightChain types
                data[0] != 0 && data[2] != 0
            },
            _ => unreachable!(),
        }
	}

	pub fn name(&self) -> String{
		match self{
			SequenceModel::IGH => "IGH-VDJ".to_string(),
			SequenceModel::IGL => "IGL-VDJ".to_string(),
			SequenceModel::IGK => "IGK-VDJ".to_string(),
			SequenceModel::TRA => "TRA-VDJ".to_string(),
			SequenceModel::TRB => "TRB-VDJ".to_string(),
			SequenceModel::TRG => "TRG-VDJ".to_string(),
			SequenceModel::TRD => "TRD-VDJ".to_string(),
			_ => unreachable!(),
		}
	}

	pub fn starts_at(&self, chain:&Chain, data:&Vec<usize> ) -> usize{
		match self{
			SequenceModel::IGH | SequenceModel::TRB | SequenceModel::TRD => {
				chain.starts_at( "HeavyChain", data) 
			},
			SequenceModel::IGL | SequenceModel::IGK | SequenceModel::TRA | SequenceModel::TRG => {
				chain.starts_at( "LightChain", data)
			},
			_ => unreachable!(),
		}
	}

}

#[derive(Clone)]
pub struct HMMcollector{
	pub states:Vec<usize>
}


impl Default for HMMcollector{
	fn default() -> Self {
        // Create a Vec<usize> with 5 zeros
        let states = vec![0; 5];
        Self { states }
    }
}

#[derive(Clone)]
pub struct HMMmodel {
	pub name:SequenceModel,
	pub collector:Vec<HMMcollector>,
}

impl HMMmodel{
	pub fn new( name:SequenceModel, size:usize) -> Self{
		let collector = vec![HMMcollector::default(); size];
		Self{
			name,
			collector,
		}
	}

	/*pub fn add (&mut self, pos:usize, value: u8) {
		if self.collector.len() < pos {
			panic!("Library was not initialized correctly - len {} is smaller than pos {}", self.collector.len(), pos );
		}
		self.collector[pos].states[HMM::char2pos(value)] +=1;
	}*/

	pub fn consume(&mut self, model: SequenceModel, start_at:usize, seq:&[u8] ) -> bool{
		if model != self.name {
			return false
		}else {
			if self.collector.len() < start_at + seq.len() {
				panic!("Library was not initialized correctly - len {} is smaller than pos {}", self.collector.len(), start_at + seq.len() );
			}
			for (pos, value) in seq.iter().enumerate(){
				for id in HMM::iupac_char2pos(*value){
					self.collector[pos+start_at].states[ id ] +=1;
				}
				
			}
			return true
		}
	}


}


pub struct VDJmodeler{

}

impl VDJmodeler{

pub fn identify_model_type(record: &str) -> Option<(SequenceModel, Chain)> {
    let name = record.to_string();

    if name.contains("IGHV") {
        return Some((SequenceModel::IGH, Chain::V));
    } else if let Some(ighd_pos) = name.find("IGHD") {
    	// Check if the 5th character after "IGHD" exists and is a digit
    	name.chars()
    		.nth(ighd_pos + 4)
    		.filter(|&c| c.is_digit(10)) // will return None if not a digit
            .map(|_| (SequenceModel::IGH, Chain::D)); {
            return Some((SequenceModel::IGH, Chain::D));
        }
        // if no - the function will return None later on.
    } else if name.contains("IGHJ") {
        return Some((SequenceModel::IGH, Chain::J));
    } else if name.contains("IGLV") {
        return Some((SequenceModel::IGL, Chain::V));
    } else if name.contains("IGLD") {
        return Some((SequenceModel::IGL, Chain::D));
    } else if name.contains("IGLJ") {
        return Some((SequenceModel::IGL, Chain::J));
    } else if name.contains("IGKV") {
        return Some((SequenceModel::IGK, Chain::V));
    } else if name.contains("IGKD") {
        return Some((SequenceModel::IGK, Chain::D));
    } else if name.contains("IGKJ") {
        return Some((SequenceModel::IGK, Chain::J));
    } else if name.contains("TRAV") {
        return Some((SequenceModel::TRA, Chain::V));
    } else if name.contains("TRAD") {
        return Some((SequenceModel::TRA, Chain::D));
    } else if name.contains("TRAJ") {
        return Some((SequenceModel::TRA, Chain::J));
    } else if name.contains("TRBV") {
        return Some((SequenceModel::TRB, Chain::V));
    } else if name.contains("TRBD") {
        return Some((SequenceModel::TRB, Chain::D));
    } else if name.contains("TRBJ") {
        return Some((SequenceModel::TRB, Chain::J));
    } else if name.contains("TRGV") {
        return Some((SequenceModel::TRG, Chain::V));
    } else if name.contains("TRGD") {
        return Some((SequenceModel::TRG, Chain::D));
    } else if name.contains("TRGJ") {
        return Some((SequenceModel::TRG, Chain::J));
    } else if name.contains("TRDV") {
        return Some((SequenceModel::TRD, Chain::V));
    } else if name.contains("TRDD") {
        return Some((SequenceModel::TRD, Chain::D));
    } else if name.contains("TRDJ") {
        return Some((SequenceModel::TRD, Chain::J));
    }

    None // If no match found, return None
}

	pub fn build_models(fasta: String) -> HMM {

		let mut reader = match parse_fastx_file(&fasta) {
        	Ok(reader) => reader,
        	Err(err) => {
            	panic!("File {fasta} Read Error: {}", err);
        	}
   		};
	    let mut models: Vec<Option<HMMmodel>> = vec![None; SequenceModel::length()];
	    let mut sequences = Vec::new();

	    // Initialize the full_matrix with the correct size
	    let mut full_matrix = vec![vec![0; 3]; SequenceModel::length()];

	    // Collect sequences and update the full_matrix with max lengths
	    while let Some(record) = reader.next() {
	    	let read = match record{
	                Ok( res ) => res,
	                Err(err) => {
	                    eprintln!("could not read from fasta:\n{err}");
	                    continue
	                }
	            };
	        let seq = read.seq().into_owned();
	        let tmp = &read.id().to_owned();
	        let acc = String::from_utf8_lossy(tmp);
	        if let Some(ids) = VDJmodeler::identify_model_type(&acc) {
	            let (model_id, chain_id) = ids;
	            full_matrix[model_id.id()][chain_id.id()] = full_matrix[model_id.id()][chain_id.id()].max(seq.len());
	            sequences.push(((model_id, chain_id), seq ));
	        }
	    }

	    // Identify SequenceModels that have sufficient data
	    let mut with_data: HashSet<SequenceModel> = HashSet::new();
	    for (id, data) in full_matrix.iter().enumerate() {
	        if let Some(seq_mod) = SequenceModel::from_index(id) {
	            if seq_mod.has_data(data) {
	                with_data.insert(seq_mod.clone());
	                let mut hmm_model = HMMmodel::new(seq_mod.clone(), data.iter().sum());
	                models[seq_mod.id()] = Some(hmm_model);
	            }
	        }
	    }
	    println!("We have found these sequences that can be modeled: {:?}", with_data);

	    // Populate the HMM models with sequence data
	    for ((model, chain), seq) in &sequences {
	        if with_data.contains(model) {
	            if let Some(hmm_model) = models[model.id()].as_mut() {
	                hmm_model.consume(
	                    model.clone(),
	                    model.starts_at(chain, &full_matrix[model.id()]),
	                    seq,
	                );
	            }
	        }
	    }

	    // Collect and return the HMM models
	    let good_models: Vec<HMMmodel> = models.into_iter().filter_map(|model| model).collect();
	    HMM::from_sequence_models(good_models)
	}
}