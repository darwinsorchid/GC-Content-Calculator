# ===================================================== DEPENDENCIES =============================================================
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from matplotlib import pyplot as plt
plt.style.use('ggplot')


# ===================================================== FUNCTIONS =================================================================
def calculate_gc(sequence):
    print(gc_fraction(sequence))

def calculate_window_gc(sequence, window, step):
    gc_parts = []
    window_count = 0
    window_list = []
    for i in range(0, len(sequence) - window +1, step):
        gc_part = gc_fraction(sequence[i: i+window])
        gc_parts.append(gc_part)
        window_count +=1
        window_list.append(f"Window {window_count}")
    return window_list, gc_parts

def create_window_plot():
    x_axis, y_axis = calculate_window_gc(sequence, window_size, step_size)
    plt.plot(x_axis, y_axis, color = "k", linestyle = "--", marker = "o")
    plt.title("Sliding Window GC Content Calculation")
    plt.xlabel("Sequence windows")
    plt.ylabel("GC content")
    plt.tight_layout()
    plt.show()

# ====================================================== LOGIC =====================================================================
sequence = Seq(input('Please enter a DNA/RNA sequence. '))

choice = input("Would you like to perform a sliding window GC content calculation? ")
if choice.lower()[0] == 'n':
    calculate_gc(sequence)
else:
    window_size = int(input("Please enter window size. "))
    step_size = int(input("Please enter step size. "))
    create_window_plot()






