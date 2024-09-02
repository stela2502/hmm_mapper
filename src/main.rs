
// main.rs
use hmm_mapper::HMM::{HMM}; // Import from lib.rs
use hmm_mapper::VDJmodeler::VDJmodeler;
use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord;

use std::io::BufReader;

use clap::Parser;


use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;


struct Seqrec{
    pub id:Vec<u8>,
    pub seq:Vec<u8>,
}

impl Seqrec{
    pub fn new (id:&[u8], seq:&[u8] ) -> Self{
        Self{
            id: id.to_vec(),
            seq: seq.to_vec(),
        }
    }

    pub fn id(&self) -> &[u8] {
        &self.id
    }
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }
}


#[derive(Parser)]
#[clap(version = "1.1.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the fasta formated IMGT database
    #[clap(short, long)]
    database: String,
    /// the fastq file you want to check for VDJ recombination events
    #[clap(short, long)]
    fastq: String,
    /// the fasta formated outfile with likely VDJ recombination evens
    #[clap(short, long)]
    outfile: String,
}

fn main() {
    let opts: Opts = Opts::parse();

    let hmm = VDJmodeler::build_models(opts.database);

    println!("Initialized HMM with {} states.", hmm.states().len());

    let fasta_path = &opts.outfile;
    let fasta_file = File::create(fasta_path).expect("Unable to create FASTA file");
    let mut fasta_writer = BufWriter::new(fasta_file);

    // Define the chunk size
    let chunk_size = 10000; // Adjust as needed

    let mut reader = match parse_fastx_file(&opts.fastq) {
        Ok(reader) => reader,
        Err(err) => {
            panic!("File {} Read Error: {}", &opts.fastq, err);
        }
    };

    let mut batch = Vec::with_capacity(chunk_size);
    let mut record_count = 0;

    while let Some(record) = &reader.next() {
        match record {
            Ok(res) => {
                batch.push( Seqrec::new( &res.id(), &res.seq()));
                record_count += 1;

                if batch.len() >= chunk_size {
                    // Process the current batch
                    process_batch(&batch, &hmm, &mut fasta_writer);
                    batch.clear(); // Clear the batch for the next set of records
                }
            }
            Err(err) => {
                eprintln!("Error reading record: {}", err);
            }
        }
    }

    // Process any remaining records in the batch
    if !batch.is_empty() {
        process_batch(&batch, &hmm, &mut fasta_writer);
    }

    println!("Processing completed. Results written to {}", fasta_path);
}

fn process_batch(batch: &[Seqrec], hmm: &HMM, fasta_writer: &mut BufWriter<File>) {
    let chunk_size = 100;
    let results : Vec<Vec<String>> = batch
    .par_chunks(chunk_size) // Specify the chunk size, e.g., 100 or another appropriate value
    .map(|chunk| {
        let mut result = Vec::<String>::with_capacity(chunk_size);
        for record in chunk.iter() {
            let seq = record.seq();
            if let Some(hmm_result) = hmm.forward_algorithm(&seq) {

                let id = String::from_utf8_lossy(&record.id()).to_string() + &format!("{:?}", hmm_result);
                let seq_str = String::from_utf8_lossy(&seq);
                let fasta_entry = format!(">{}\n{}\n", id, seq_str);
                result.push( fasta_entry );
            }
        }
        return result
    })
    .collect();

    let mut it = 0;
    for result_vec in results {
        it += result_vec.len();
        for result in &result_vec{
            //println!( "{:?}", result);
            write!(fasta_writer, "{}", result).expect("Failed to write to FASTA file");
        }
    }
    println!("Batch processed -> {it} potential VDJ reads found.");
}

/*fn main() {

    let opts: Opts = Opts::parse();

    let hmm = VDJmodeler::build_models( opts.database );

    println!("Initialized HMM with {} states.", hmm.states().len());

    let mut reader = match parse_fastx_file(&opts.fastq) {
        Ok(reader) => reader,
        Err(err) => {
            panic!("File {} Read Error: {}",&opts.fastq, err);
        }
    };

    // Iterate over the records
    while let Some(record) = reader.next() {
        let read = match record{
            Ok( res ) => {
                if let Some( hmm_result ) = hmm.forward_algorithm(&res.seq()){
                    //if hmm_result.iter().take(hmm_result.len() - 1).any(|x| x.1 < 1e-08){
                        /*println!( "This {:?} is what the HMM tells me for \n{:?}\n{:?}",
                            hmm_result,
                            String::from_utf8_lossy( &res.id() ),
                            String::from_utf8_lossy( &res.seq() ),
                        );*/
                        println!( ">{}\n{}", 
                            String::from_utf8_lossy( &res.id() ),
                            String::from_utf8_lossy( &res.seq() ),
                            );
                    //}
                    
                }/*else {
                    eprintln!("Couln not read this sequence: {}",String::from_utf8_lossy( &res.seq() ) )
                }*/
                
            },
            Err(err) => {
                eprintln!("could not read from fasta:\n{err}");
                continue
            }
        };
    }

}

*/