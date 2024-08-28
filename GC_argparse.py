# ===================================================== DEPENDENCIES ==============================================================
import argparse
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from matplotlib import pyplot as plt
plt.style.use('seaborn-v0_8')

# ===================================================== ARGPARSE SETUP =============================================================
# Create parser instance
parser = argparse.ArgumentParser(description='Calculate GC content of sequence')
# Add positional arguments
parser.add_argument('seq', nargs= '+', help = 'Enter a single or multiple DNA/RNA sequences of the same length')

# Parse command-line arguments that were passed to the script and store them as attributes to args object
args = parser.parse_args()

# Convert sequences from strings to Seq objects
sequence_list= [Seq(seq) for seq in args.seq]

# ===================================================== FUNCTIONS ==================================================================
def calc_gc(sequence_list):
    ''' Calculates GC content of sequences and returns list of GC content of each sequence.
    '''
    gc_list = []
    for i in sequence_list:
        gc_list.append(round(gc_fraction(i), 2))
    return gc_list


def create_GC_barplot(sequence_list, gc_list):
    '''Function that creates barplot of the GC value calculated (not with SWAN) of multiple sequences.
    '''
    # Convert sequences in sequence_list from Seq objects to strings
    x_values = [str(seq) for seq in sequence_list]

    fig, ax = plt.subplots()
    
    # Create the bar plot
    bar_container = ax.bar(x_values, gc_list, width = 0.3, color = '#3937b3')
    
    # Set plot titles and labels
    ax.set_title("Sequence GC Content Analysis")
    ax.set_xlabel("Sequences")
    ax.set_ylabel("GC content")

    # Add labels on the bars
    ax.bar_label(bar_container)

    # Set axis limits if only one sequence is present
    if len(sequence_list) == 1:
        ax.set_xlim(-0.5, 0.5)

    # Style the spines (borders) of the plot
    ax.spines['top'].set_color('0.5')
    ax.spines['right'].set_color('0.5')
    ax.spines['bottom'].set_color('0.5')
    ax.spines['left'].set_color('0.5')
    
    # Display the plot
    plt.show()

# ===================================================== LOGIC ===============================================================
# Calculate gc_content of passed in sequences
y_values = calc_gc(sequence_list)

# Create barplot of GC content values
create_GC_barplot(sequence_list, y_values)
