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
parser.add_argument('seq', nargs= '+', help = 'Enter a single or multiple DNA/RNA sequences of the same length!')
parser.add_argument('-w', '--window', type = int, help = 'Enter the window size if you want to calculate GC content with sliding window analysis.')
parser.add_argument('-s', '--step', type = int, help = 'Enter step size if you want to calculate GC content with sliding window analysis.' )

# Parse command-line arguments that were passed to the script and store them as attributes to args object
args = parser.parse_args()

# Convert sequences from strings to Seq objects
sequence_list= [Seq(seq) for seq in args.seq]

# ===================================================== FUNCTIONS ==================================================================

# ------------------------------------- CALCULATE GC CONTENT WITHOUT SLIDING WINDOW ------------------------------------------------
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

# ---------------------------------------CALCULATE GC CONTENT WITH SLIDING WINDOW -------------------------------------------
def calculate_window_gc(seq, window, step):
    '''Function that calculates the GC content of sequenceS with the sliding window method.
    Returns a tuple of a list of the windows the sequence has been split to and a list with 
    the respective GC values of every window analyzed.

    [Sliding window analysis (SWAN) is a commonly used method for studying the properties of molecular sequences:
    data are plotted as moving averages of a particular criterion, such as the number of nucleotide changes,
    for a window of a certain length slid along a sequence or sequence alignment (Tajima, 1991).
    SWAN developed to reveal patterns of variability along nucleotide sequences.]
    '''
    gc_parts = []
    window_count = 0
    window_list = []

    # iterate through windows in a single sequence and calculate GC content for each window
    for i in range(0, len(seq) - window +1, step):
        gc_part = gc_fraction(seq[i: i+window])
        gc_parts.append(gc_part)
        window_count +=1
        window_list.append(f"Window {window_count}")
    return window_list, gc_parts


def multi_window_plot(x_windows, seq_dict):
    '''Function that plots the SWAN calculated GC content of multiple sequences passed in by the user. Window and step size are the same
    for all the sequences analyzed and they are defined by the user. The function takes in as an argument a dictionary with the sequences 
    as keys and a list of GC values of each window as values.
    '''
    colors_list = ['b', 'g', 'm','r', 'c', 'y', 'k']
    for num, (key, value) in enumerate(seq_dict.items()):
        color = colors_list[num % len(colors_list)]
        plt.plot(x_windows[0], value, color = color, marker = "o", linestyle = "--", label = f"{key}")
    plt.title("Sliding Window GC Content Calculation")
    plt.xlabel("Sequence Windows")
    plt.ylabel("GC content")
    plt.legend()
    plt.tight_layout()
    plt.show()

# ===================================================== LOGIC ===============================================================

# Calculate GC content of one or more RNA/DNA sequences with sliding window analysis
if args.window is not None and args.step is not None:
    windows = []
    gc_content_list = []
    
    # Iterate through each given sequence
    for seq in sequence_list:
        window_list, gc_content = calculate_window_gc(seq, args.window, args.step)
        windows.append(window_list)
        gc_content_list.append(gc_content)
    
    # Create dictionary for given sequences and number of respective windows
    pairs = zip(sequence_list, gc_content_list)
    seq_dict = dict(pairs)
    multi_window_plot(windows, seq_dict)
elif args.window is not None and args.step is None:
    print("Please provide step size")
elif args.window is None and args.step is not None:
    print("Please provide window size")

# Calculate GC content of one or more RNA/DNA sequences without sliding window analysis
else:
    # Calculate gc_content of passed in sequences
    y_values = calc_gc(sequence_list)

    # Create barplot of GC content values
    create_GC_barplot(sequence_list, y_values)