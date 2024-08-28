//VDHmodeler.rs

use crate::HMM::HMM;
use bio::io::fastq;


pub enum LightChain{
	V,
	J
}


pub enum HeavyChain{
	V,
	D,
	J
}

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


	pub fn have_data(&self, data:&Vec<usize> ) -> bool{
		match self {
            SequenceModel::IGH(_) | SequenceModel::TRA(_) | SequenceModel::TRG(_) => {
                // Check for HeavyChain types
                data[0] != 0 && data[1] != 0 && data[2] != 0
            },
            SequenceModel::IGL(_) | SequenceModel::IGK(_) | SequenceModel::TRB(_) | SequenceModel::TRD(_) => {
                // Check for LightChain types
                data[0] != 0 && data[1] != 0
            },
        }
	}

}

#[derive(Clone)]
struct HMMcollector{
	states:Vec<usize>
}


impl Default for HMMcollector{
	fn default() -> Self {
        // Create a Vec<usize> with 5 zeros
        let states = vec![0; 5];
        Self { states }
    }
}


pub struct HMMmodel {
	pub name:SequenceModel,
	collector:Vec<HMMcollector>,
}

impl HMMmodel{
	pub fn new( name:SequenceModel, size:usize) -> Self{
		let collector = vec![HMMcollector::default(); size];
		Self{
			name,
			collector,
		}
	}

	pub fn add (&mut self, pos:usize, value: char) {
		if self.collector.len() < pos {
			panic!("Library was not initialized correctly - len {} is smaller than pos {}", self.collector.len(), pos );
		}
		self.collector[pos].states[HMM::char2pos(value)] +=1;
	}


}



pub fn build_models( fasta : String ) -> Vec<HMMmodel> {

	let reader = fastq::Reader::from_file(file_path).unwrap();
    let mut sequences = Vec::new();

    let full_matrix = vec[vec[0;3], SequenceModel.length()];

    for record in reader.records() {
        let record = record.unwrap();
        let seq = record.seq().to_owned();
        if Some((ids)) = identify_model_type(&record.id().to_string()) {
        	full_matrix[ ids.0.id()][ids.1.id() ] = full_matrix[ ids.0.id()][ids.1.id() ].max(seq.len());
        	sequences.push( ( ids, seq) );
        }
    }
    // check if we have sequences for all necessary parts for each of the SequenceModel(s)
    let mut with_data: Vec<HMMmodel>;
    for id in 1:full_matrix.len(){
    	if SequenceModel::from_index(id).has_data( &full_matrix[id] ){
    		
    	}
    }

    sequences
}