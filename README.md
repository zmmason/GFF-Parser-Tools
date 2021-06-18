# GFF-Parser-Tools (Version 3.2)
Bioinformatic Analysis Tools (transcriptomics, genomics, proteomics), File Conversion, Gene Expression Profiles, FASTA Generator, Gene-Ontology, SNPs, CDSs/Introns/Transcripts, GO-terms, Pfam parser, FASTA reader, Molar Extinction, Charge, NT-AA Composition, Isoelectric Point,  Mass Extinction, Molecular Weight, Codon Composition,  ORF Reader,  Reverse Compliment.

# GFF_Parser.py
This program is designed to take an input .txt file containing significant gene/feature labels and search a reference GFF file for their corresponding parent sequence ID and region index. The genes are then extracted from the provided FASTA file using the parent sequence titles and indexes associated with each gene. The single-line FASTA formatted output file is designed to be used in database searches. NOTE: There are file size caps on many databases such as BLAST. 

Input (commandline): 
- Feature List (.txt)
- Refference GFF file (.GFF/.GFF3)
- Asscociated refference FASTA file (.fa/.fas/.txt)

Output: 
- nt-base FASTA (.fa).

# Sequence_Analysis.py
This program is a library of tools used to analyze genomic sequences:
  - FASTA Reader
  - Molar Extinction
  - Charge
  - NT-AA Composition 
  - Theoretical Isoelectric Point  
  - Mass Extinction 
  - Molecular Weight
  - Codon Composition
  - ORF Reader
  - Reverse Compliment Generator

# Contributors
- Zachary M Mason (zmmason@ucsc.edu) - GFF_Parser.py, Sequence_Analysis.py
- Justin Chan (jumchan@ucsc.edu) - Sequence_Analysis.py
