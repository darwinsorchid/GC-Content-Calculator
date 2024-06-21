# ===================================================== DEPENDENCIES =============================================================
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from matplotlib import pyplot as plt
plt.style.use('ggplot')


# ===================================================== FUNCTIONS =================================================================
def calculate_gc(sequence):
    ''' Simple function that calculates GC content of a single sequence and rounds the floating point number to the third decimal place.
    '''
    print(round(gc_fraction(sequence), 3))


def multiple_gc(seq_choice):
    ''' Function that takes in multiple sequences as input from user and calculates GC content. Returns tuple of a list of the sequences
    and a list of the GC contents of the sequences.
    '''
    sequences = []
    gc_list = []
    print("All sequences must be of the same length.")
    for i in range(seq_choice):
        sequence = Seq(input("Please enter a DNA/RNA sequence. "))
        sequences.append(f"Seq {i+1}")
        gc_list.append(round(gc_fraction(sequence), 3))
    return sequences , gc_list


def calculate_window_gc(sequence, window, step):
    '''Function that calculates the GC content of a single sequence with the sliding window method. The window and step size are passed in
    as positional arguments and the user defines their values accordingly. Returns a tuple of a list of the windows the sequence has been 
    split to and a list with the respective gc values of every window analyzed.

        [Sliding window analysis (SWAN) is a commonly used method for studying the properties of molecular sequences:
        data are plotted as moving averages of a particular criterion, such as the number of nucleotide changes,
        for a window of a certain length slid along a sequence or sequence alignment (Tajima, 1991).
        SWAN developed to reveal patterns of variability along nucleotide sequences.]
    '''
    gc_parts = []
    window_count = 0
    window_list = []
    for i in range(0, len(sequence) - window +1, step):
        gc_part = gc_fraction(sequence[i: i+window])
        gc_parts.append(gc_part)
        window_count +=1
        window_list.append(f"Window {window_count}")
    return window_list, gc_parts


def create_window_plot(sequence):
    '''Function that takes plots the SWAN calculated content of a single sequence passed in by the user.
    '''
    x_axis, y_axis = calculate_window_gc(sequence, window_size, step_size)
    plt.plot(x_axis, y_axis, color = "k", linestyle = "--", marker = "o", label = f"{sequence}")
    plt.title("Sliding Window GC Content Calculation")
    plt.xlabel("Sequence Windows")
    plt.ylabel("GC content")
    plt.legend()
    plt.tight_layout()
    plt.show()


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


def create_GC_barplot(x, y):
    '''Function that creates barplot of the GC value calculated (not with SWAN) of multiple sequences.
    '''
    fig, ax = plt.subplots()
    bar_container = ax.bar(x, y, width = 0.4, color = '#C04000')
    ax.set_title("Sequence GC Content Analysis")
    ax.set_xlabel("Sequences")
    ax.set_ylabel("GC content")
    ax.bar_label(bar_container)
    ax.spines['top'].set_color('0.5')
    ax.spines['right'].set_color('0.5')
    ax.spines['bottom'].set_color('0.5')
    ax.spines['left'].set_color('0.5')
    fig.tight_layout()
    plt.show()

    
# ====================================================== LOGIC =====================================================================
print("Welcome to GC calculator!")
seq_choice = int(input("Please enter the number of the sequences you want to analyze today. "))
if seq_choice == 1:
    sequence = Seq(input('Please enter the DNA/RNA sequence. '))
    choice = input("Would you like to perform a sliding window GC content calculation? ")
    if choice.lower()[0] == 'n':
        calculate_gc(sequence)
    else:
        window_size = int(input("Please enter window size. "))
        step_size = int(input("Please enter step size. "))
        create_window_plot(sequence)
        
else:
    mult_choice = input("Would you like to perform a sliding window GC content calculation?")
    if mult_choice.lower()[0] == 'n':
        x_values , y_values = multiple_gc(seq_choice)
        print(x_values)
        print(y_values)
        create_GC_barplot(x_values, y_values)
        
    else:
        print("Note! All sequences will be analyzed using the same window and step size.")
        window_size = int(input("Please enter window size. "))
        step_size = int(input("Please enter step size. "))
        sequence_list = []
        windows = []
        gc_content_list = []
        print("All sequences must be of the same length.")
        for i in range(seq_choice):
            new_seq = Seq(input("Please enter a DNA/RNA sequence. "))
            sequence_list.append(new_seq)
        for seq in sequence_list:
            window , gc_content = calculate_window_gc(seq, window_size, step_size)
            windows.append(window)
            gc_content_list.append(gc_content)
        # create dictionary of sequences and gc_values
        pairs = zip(sequence_list, gc_content_list)
        seq_dict = dict(pairs)
        multi_window_plot(windows, seq_dict)


        




