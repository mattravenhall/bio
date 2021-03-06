### General Bioinformatics Module
This is a general bioinformatics Python module in progress. Whilst all core functions should work as intended, they are not neccessarily optimised. Redevelopment to resolve this is ongoing but irregular. Whilst bio generally supports both Python 2 & 3, Python 3 should be used to ensure compatibility.

#### biocore.py
Completed miscellaneous functions, most non-specialised functions should be found here.

Usage: biocore [command] [arguments]

Commands:

  - predictMT: Rough melting temperature prediction.
  - findkmers: Find kmers of given length within a fasta.
  - translate: Translate from DNA/RNA to Protein, auto-detects.
  - transcribe: Transcribe from RNA/DNA to DNA/RNA, auto-detects.
  - complement: Find the complement of a DNA sequence.
  - calcHamming: Determine the Hamming distance between two sequences.
  - simCleave: Simulate cleavage of a sequence by a given enzyme.
  - simCleaveMulti: simCleave for multiple sequences provided as a fasta/q.
  - simPCR: Predict PCR fragments of a given sequence and two primers.
  - simPCRMulti: simPCR for multiple sequences provided as a fasta/q.
  - scaffToContigs: Convert single scaffold genome to contigs.
  - findMotif: Given a motif, find start positions in fasta file or sequence.
  - AAchange: Given a SNP and gene sequence, predict amino acid change.
  - BPtoAA: Convert a genomic position to an amino acid position.
  - AAtoBP: Convert an amino acid position to genomic positions.
  - consensus: Find a consensus sequence for a multi-contig fasta/q.
  - profile: Produce a profile matrix for a given multi-contig fasta/q.
  - contigExtract: Write out specified contigs from a fasta/q file.
  - seqExtract: Write out a specific sequence from a fasta/q file.
  - fasta: Create an example fasta sequence for testing.
  - palindrome: Given a sequence, determine the longest palindrome present.

NB: This list is likely incomplete, further internal functions may exist and can be called within a Python intrepreter.

#### biostats.py
Functions of statistical importance, this script is primarily intended to return overview statistics for newly assembled genomes but may later be adapted for more general evaluation.

Usage: biostats [command] [arguments]

Commands:
  - full: Overview statistics for a given fasta.
  - topGC: Returns the contig/strain with the highest GC.

### Additional .py files
Any functions additional to those listed above should be considered as under active development, they may therefore produce unexpected/unintended results. Use with additional caution.

#### Acknowledgements
Thanks to Álvaro Abella Bascarán for style corrections.
