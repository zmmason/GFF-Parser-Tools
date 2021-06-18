# GFF-Parser-Tools (Version 3.2)

Pull target features from your reference FASTA using the associated GFF/GFF3 to create a single-line FASTA file containing your target features.

Title format: feature-title|parent-seq-name|NT|feature|start:stop

Bioinformatic Analysis Tools (transcriptomics, genomics, proteomics), File Conversion, Gene Expression Profiles, FASTA Generator, Gene-Ontology, SNPs, CDSs/Introns/Transcripts, GO-terms, Pfam parser, FASTA reader, Molar Extinction, Charge, NT-AA Composition, Isoelectric Point,  Mass Extinction, Molecular Weight, Codon Composition,  ORF Reader,  Reverse Compliment.

# GFF_Parser.py
This program is designed to take an input .txt file containing significant gene/feature labels and search a reference GFF file for their corresponding parent sequence ID, reading frame, stand, and region index. The features are then extracted from the provided FASTA file using the parent sequence titles and indexes associated with each gene from the appropraite stand & frame. The single-line FASTA formatted output file is designed to be used in database searches/comparisons. 
NOTE: Only handles GFF entries with 1 transcript ('.t1'). Additional transcripts will be ignored unless the target feature is a 'gene'.' 

Input (terminal interaction): 
- Feature List (.txt)
- Refference GFF file (.GFF/.GFF3)
- Asscociated refference FASTA file (.fa/.fas/.txt)
- Target Feature (gene/intron/exon/CDS/transcript/match/rRNA/tRNA/mRNA)
- Output file name

Output: 
- nt-base FASTA

# Sequence_Analysis.py
This program is a library of tools used to analyze genomic sequences:
  - FASTA Reader
  - Molar Extinction
  - Charge
  - NT/AA Composition 
  - Theoretical Isoelectric Point  
  - Mass Extinction 
  - Molecular Weight
  - Codon Composition
  - ORF Reader
  - Reverse Compliment Generator

# Contributors
- Zachary M Mason (zmmason@ucsc.edu) - GFF_Parser.py, Sequence_Analysis.py
- Justin Chan (jumchan@ucsc.edu) - Sequence_Analysis.py
