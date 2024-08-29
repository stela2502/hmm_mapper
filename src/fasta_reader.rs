use std::fs::File;
use std::io::{BufRead, BufReader, Read, Result};
use flate2::read::GzDecoder;
use core::fmt;
use std::error::Error;

/// Structure to hold a single Fasta record
#[derive(Debug)]
pub struct FastaRecord {
    pub header: String,
    pub sequence: String,
}

impl FastaRecord {
    /// Create a new FastaRecord
    pub fn new(header: String, sequence: String) -> Self {
        FastaRecord { header, sequence }
    }

    pub fn seq(&self) -> &str {
        &self.sequence
    }

    pub fn id(&self) -> &str {
        &self.header
    }
}

impl fmt::Display for FastaRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Write the header
        writeln!(f, "{}", self.header)?;

        // Write the sequence in 60 bp lines
        for chunk in self.sequence.as_bytes().chunks(60) {
            writeln!(f, "{}", std::str::from_utf8(chunk).unwrap())?;
        }

        Ok(())
    }
}

/// Structure to handle reading Fasta files
pub struct FastaReader {
    records: Vec<FastaRecord>,
}

impl FastaReader {
    /// Open a Fasta file, read its content, and create a FastaReader
    pub fn from_file(filename: &str) -> Result<Self> {
        let file = File::open(filename)?;
        let mut buf_reader = BufReader::new(file);

        // Check if the file is gzipped
        let mut reader: Box<dyn BufRead> = if Self::is_gzipped(&mut buf_reader) {
            Box::new(BufReader::new(GzDecoder::new(buf_reader)))
        } else {
            Box::new(buf_reader)
        };

        let mut content = String::new();
        reader.read_to_string(&mut content)?;

        let records = Self::parse_records(&content);
        Ok(FastaReader { records })
    }

    /// Check if a file is gzipped by reading its magic number
    fn is_gzipped<R: Read>(reader: &mut R) -> bool {
        let mut magic_number = [0; 2];
        match reader.read_exact(&mut magic_number) {
            Ok(_) => magic_number == [0x1f, 0x8b],
            Err(_) => false,
        }
    }

    /// Parse the FASTA file content into records
    fn parse_records(content: &str) -> Vec<FastaRecord> {
        let mut records = Vec::new();
        let mut lines = content.lines();

        let mut sequence: String = "".to_string();
        let mut header: Option<String> = None;

        while let Some(data) = lines.next() {

            if data.starts_with('>') {
                if let Some(id) = header{
                    records.push ( FastaRecord::new( id, sequence ) );
                }
                header = Some(data.to_string());
                sequence = "".to_string();
            }else {
                sequence += &data
            }
        }

        if let Some(id) = header{
            records.push ( FastaRecord::new( id, sequence ) );
        }

        records
    }
}

impl Iterator for FastaReader {
    type Item = FastaRecord;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.pop()
    }
}