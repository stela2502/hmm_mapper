# hmm_mapper

**hmm_mapper** is a tool designed to detect VDJ recombination events in DNA sequences using Hidden Markov Models (HMMs). While it is a work in progress, it aims to be a valuable component for mapping and analyzing VDJ sequences.

## What It Does

**hmm_mapper** performs the following tasks:
1. **VDJ Detection Module**: Utilizes information from the IMGT database to gather all possible VDJ recombination events.
2. **HMM Scoring**: Constructs an HMM based on the collected VDJ data and scores sequences for potential VDJ recombination.
3. **Sequence Analysis**: Accepts sequences in Fastq or Fasta format and returns those likely to contain VDJ recombination events, along with their accession and sequence data.

While currently in its early stages, **hmm_mapper** has the potential to become a powerful tool when integrated into a larger mapping framework.

## Installation

### Prerequisites

- Ensure you have the [Rust compiler](https://www.rust-lang.org/tools/install) installed.

### Installation Steps

Clone the repository and build the project:
```bash
git clone https://github.com/stela2502/hmm_mapper
cd hmm_mapper
cargo build --release
```

Alternatively, you can install it directly from the GitHub repository:
```bash
cargo install --git https://github.com/stela2502/hmm_aligner
```

## Usage

To see the available options and usage instructions, run:
```bash
hmm_mapper -h
```

### Example Command

```bash
hmm_mapper --database <DATABASE> --fastq <FASTQ>
```

### Options

- `-d, --database <DATABASE>`: Path to the IMGT database in Fasta format.
- `-f, --fastq <FASTQ>`: Path to the Fastq file you want to analyze for VDJ recombination events.
- `-h, --help`: Displays help information.
- `-V, --version`: Displays version information.

## Work in Progress

Please note that **hmm_mapper** is still under development, and its full capabilities and usability are yet to be fully assessed.





