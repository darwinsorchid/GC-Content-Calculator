# GC-Content Calculator 

## Description 
* Tool for calculating and plotting the GC content of one or multiple sequences taken in as user input. 
* User can choose between calculating total GC content of each sequence or performing sliding window analysis (SWAN) for the calculation of the GC content for each window,
  to reveal patterns of variability along nucleotide sequences.

## Sequence input options
- [GC_user_input] DNA/RNA sequences, as well as window and step size in case of sliding window analysis, are passed in via user input.
- [GC_argparse] DNA/RNA sequences, as well as window and step size in case of sliding window analysis, are passed in via the command line: At least one sequence is required as a positional argument,
whereas window and step size are optional arguments.
- [GC_argparse_file] DNA/RNA sequences are parsed from different file formats passed in via the command line as optional arguments. The user can choose between the following file formats:
  1. FASTA (.fasta)      
  2. FASTA nucleid acid (.fna) 
  3. GenBank (.gb)
  4. Aligned FASTA (.fa)
  5. ClustalW (.aln)
  
  _**NOTE** : Sliding window analysis option is not available for calculating GC content of sequences from files._

## Expected outputs
[GC_user_input]
* Single sequence / Simple calculation    : Floating point number of GC value
* Single sequence / SWAN                  : Simple line plot of GC values
* Multiple sequences / Simple calculation : Barplot of GC values
* Multiple sequences / SWAN               : Line plot with multiple y values

_**NOTE** : Sequences in figures are represented as they are passed in._

[GC_argparse]
* Single sequence / Simple calculation    : Barplot of GC value
* Single sequence / SWAN                  : Simple line plot of GC values
* Multiple sequences / Simple calculation : Barplot of GC values
* Multiple sequences / SWAN               : Line plot with multiple y values
  
_**NOTE** : Sequences in figures are represented as they are passed in._

[GC_argparse_file]
* Single sequence    : Barplot of GC value
* Multiple sequences : Barplot of GC values
  
_**NOTE** : Sequences in figures are represented as their respective accession numbers / species names depending on the file format._

## Libraries
| [GC_user_input] | [GC_argparse] | [GC_argparse_file] |
| --------------- |:-------------:| ------------------:|
|   matplotlib    |    argparse   |       argparse     |
|   biopython     |   matplotlib  |      matplotlib    |
|                 |   biopython   |      biopython     |
|                 |               |         sys        |
