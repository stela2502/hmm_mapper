
// main.rs
use hmm_mapper::HMM::{HMM}; // Import from lib.rs
use hmm_mapper::VDJmodeler::VDJmodeler;
use needletail::parse_fastx_file;

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

