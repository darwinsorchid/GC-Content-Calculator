# ===================================================== DEPENDENCIES ==============================================================
import argparse
import sys
from Bio import SeqIO
from Bio import Align
from Bio.SeqUtils import gc_fraction
from matplotlib import pyplot as plt
plt.style.use('seaborn-v0_8')

# ===================================================== ARGPARSE SETUP =============================================================
# Create parser instance
parser = argparse.ArgumentParser(description = 'Calculate GC content of DNA/RNA sequences from files of different formats.')

# ----------------------------- Single sequence file formats as optional arguments ------------------------------
# FASTA nucleic acid files
parser.add_argument('-fna', nargs='+', type = argparse.FileType('r', encoding = 'UTF-8'), help = 'Provide .fna file')

# FASTA files
parser.add_argument('-fasta', nargs = '+', type = argparse.FileType('r', encoding = 'UTF-8'), help = 'Provide .fasta file')

# GenBank files
parser.add_argument('-gb', nargs = '+', type = argparse.FileType('r', encoding = 'UTF-8'), help = 'Provide .gb file')


# ------------------------------ Sequence alignment file formats as optional arguments ------------------------------
# Aligned FASTA files
parser.add_argument('-fa', type = argparse.FileType('r', encoding = 'UTF-8'), help = 'Provide .fa file')

# ClustalW files
parser.add_argument('-aln', type = argparse.FileType('r', encoding = 'UTF-8'), help = 'Provide .aln file')

# Parse command-line arguments that were passed to the script and store them as attributes to args object
args = parser.parse_args()

# Check if no files were provided
if not any([args.fna, args.fasta, args.gb, args.fa, args.aln]):
    print("Provide a file to continue")
    sys.exit(1)  # Exit the script with a status of 1 (indicating an error)

elif args.fna is not None:
    record_list =[SeqIO.read(file, "fasta") for file in args.fna]            # Make list of SeqRecord objects from file
    sequence_list = [record.seq for record in record_list]                   # Make list of sequences as Seq objects
    accession_num_list = [record.id.split('|')[3] for record in record_list] # Make list of accession numbers of sequences

elif args.fasta is not None:
    record_list =[SeqIO.read(file, "fasta") for file in args.fasta]          
    sequence_list = [record.seq for record in record_list]                   
    accession_num_list = [record.id for record in record_list]               

elif args.gb is not None:
    record_list =[SeqIO.read(file, "genbank") for file in args.gb]          
    sequence_list = [record.seq for record in record_list]
    accession_num_list = [record.id for record in record_list]

elif args.fa is not None:
    alignment = Align.read(args.fa, 'fasta')
    record_list = alignment.sequences
    sequence_list = [record.seq for record in record_list]
    accession_num_list = [record.id for record in record_list]     # Species name instead of accession number for aligned FASTA files

elif args.aln is not None:
    alignment = Align.read(args.aln, 'clustal')
    record_list = alignment.sequences
    sequence_list = [record.seq for record in record_list]
    accession_num_list = [record.id.split('|')[3] for record in record_list]

# ===================================================== FUNCTIONS ==================================================================
def calc_gc(sequence_list):
    ''' Calculates GC content of sequences and returns list of GC content of each sequence.
    '''
    gc_list = []
    for i in sequence_list:
        gc_list.append(round(gc_fraction(i), 2))
    return gc_list


def create_GC_barplot(accession_num_list, gc_list):
    '''Function that creates barplot of the GC value calculated of multiple sequences without sliding window analysis.
    '''
    fig, ax = plt.subplots()

    # Create the bar plot
    bar_container = ax.bar(accession_num_list, gc_list, width = 0.3, color = '#3937b3')
    
    # Set plot titles and labels
    ax.set_title("Sequence GC Content Analysis")
    ax.set_xlabel("Sequences")
    ax.set_ylabel("GC content")

    # Add labels on the bars
    ax.bar_label(bar_container)

    # Set axis limits if only one sequence is present
    if len(accession_num_list) == 1:
        ax.set_xlim(-0.5, 0.5)

    # Style the spines (borders) of the plot
    ax.spines['top'].set_color('0.5')
    ax.spines['right'].set_color('0.5')
    ax.spines['bottom'].set_color('0.5')
    ax.spines['left'].set_color('0.5')
    
    # Display the plot
    plt.show()


# Calculate GC content of sequences from file
gc_content_list = calc_gc(sequence_list)

# Create barplot of GC content
create_GC_barplot(accession_num_list, gc_content_list)
