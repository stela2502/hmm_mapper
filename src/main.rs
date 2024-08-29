
// main.rs
use hmm_mapper::HMM::{HMM}; // Import from lib.rs
use hmm_mapper::VDJmodeler::VDJmodeler;
use needletail::parser::FastqReader;
use std::io::BufReader;

use clap::Parser;


#[derive(Parser)]
#[clap(version = "1.1.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the fasta formated IMGT database
    #[clap(short, long)]
    database: String,
    /// the fastq file you want to check for VDJ recombination events
    #[clap(short, long)]
    fastq: String,
}


fn main() {

    let opts: Opts = Opts::parse();

    let hmm = VDJmodeler::build_models( opts.database );

    println!("Initialized HMM with {} states.", hmm.states().len());

    let reader = BufReader::new(opts.fastq);

    // Create a FastqReader from the buffered reader
    let mut fastq_reader = FastqReader::new(reader);

     // Iterate over the records
    for record in fastq_reader.records() {
        match record {
            Ok(rec) => {
                // Process the record
                let id = rec.id();
                let sequence = rec.seq();
                let quality = rec.qual();

                println!("ID: {}", id);
                println!("Sequence: {}", std::str::from_utf8(sequence).unwrap());
                println!("Quality: {}", std::str::from_utf8(quality).unwrap());
            }
            Err(e) => {
                eprintln!("Error reading FASTQ record: {}", e);
            }
        }
    }

}

