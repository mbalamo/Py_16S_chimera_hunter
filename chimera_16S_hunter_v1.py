import sys, os, re
import pandas as pd
import matplotlib.pyplot as plt

# Purpose of this script is to take a FASTA file that contains long 16S sequences, correct the line wrap to have 2-lines per 
# entry, and then search each contig to check if both ends match to the same genus or not. If not, that contig is chimeric. 
# The FASTA contigs are derived from the zymo microbial control sample, and have known 16S sequences. We can use this sample
# and this "hunter" script to determine if our process or reagents for library prep is introducing chimeric molecules, which 
# we want to avoid. It also looks for truncated contigs, which can be another symptom of a poor prep process. It also looks
# for contigs that do not match with the zymo sequences for 16S, suggesting off-target amplification or low quality reads
# with a lot of errors. One the weakness here is that the script assumes/requires perfect matchs to the 16S V1 and V9 target
# ends. But this is also a biological nescessity, because 16S across genus can be highly similiar to nearly ideintical. A
# single mismatch can be a real difference for calling one genus over another. Therefore, we cannot accept mismatches in our 
# anchor matching since mismatches can be real events that differentiate genus.
# 
# After correcting the FASTA file to have 2 lines per entry, the script first defines the anchors, V1(start) and V9(end) of 
# 16S for all 8 zymo genus. It then reads the FASTA file, one entry at a time, and checks each anchor against the seqeunce 
# looking for the start and end matches using the key:value to assign a genus(key) based on the anchor(values) that match.
# As the ends are named to a genus (or undefined if no match is found), the end names are compared to see if they are from 
# the same genus (a complete 16S) or mis-match (chimera). Truncated contigs are counted but not serched for matches. The 
# counts are summed over the file, then a summary txt is saved along with a by_count and by_percentage bar plots to show
# the final counts for all categories examined.

# highly similar V1 so pick sites unique across all 8 genus in the zymo mix.
Bacillus_V1 = "CGAGCGGACAGATGGGAGCTTGCTCCCTGAT"
Enterococcus_V1 = "GAACGCTTCTTTCCTCCCGAGTGCTTGCACTCATTTGGAAA"
Lactobacillus_V1 = "TTGGCCCAATTGATTGATGGTGCTTGCACCTGATTGATTTTGGTC"
Listeria_V1 = "CGAACGAACGGAGGAAGAGCTTGCTCTTCCAAA"
Pseudomonas_V1 = "TCGAGCGGATGAAGGGAGCTTGCTCCTGGA"
Staphylococcus_V1 = "GAGCGAACGGACGAGAAGCTTGCTTCTCTG"
# these have 2 types of 16S -- pick seqeunce match to account for both
Escherichia_V1 = "AGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGAC"
Salmonella_V1 = "AGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGAC"

# find unique ends as close to the end as possible, unique acrosss all in the zymo mix
Bacillus_V9 = "CCAGCCGCCGAAGGTGGGACAGATGATTGG"
Enterococcus_V9 = "CCAGCCGCCTAAGGTGGGATAGATGATTG"
Staphylococcus_V9 = "CTAGCCGTCGAAGGTGGGACAAATGATTGGGGT"
Lactobacillus_V9 = "TTTGTAACACCCAAAGTCGGTGGGGTAACCTTTTAGG"
Listeria_V9 = "CCAGCCGCCGAAGGTGGGACAGATAATTG"
Pseudomonas_V9 = "GACGGTTACCACGGAGTGATTCATGACTGGGG"
# these have 2 types of 16S per genus - pick seqeunce match to account for both, around 1193 bp
Escherichia_V9 = "TTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGA"
Salmonella_V9 = "CAGCGATTAGGTCGGGAACTCAAAGG"

# -------------------------------------------
def FASTA_fixwrap(filename):
# Purpose of this function is to fix newlines in a FASTA file so it becomes 2-line file.
    try:
        # Open the FA fasta file.
        input_file = open(filename,'r')
        # Create a new destination file as a temp file of the corrected FASTA file.
        newname = "fixed_temp.fasta"
        output_file = open(newname, 'w')
        # The very first line is a header and end in \n so write it directly to a temp file
        output_file.write(input_file.readline())
        # now go over the rest of the file
        for each_line in input_file:
            if (each_line[:1] == ">"):
                # is the start of a new entry, so write it direct. This line guarneteed to be one line ending in \n
                # but preface with a newline to account for the seq of a previous entry NOT having a newline element.
                output_file.write("\n" + each_line)
            else:
                # trim off the newline at the end of this line then write it out
                output_file.write(each_line[:-1])
        input_file.close()
        output_file.close()
        return newname

    except Exception as e:
        print ("An error was generated when attempting to modify the FASTA file.")
        exit ()
# -------------------------------------------

def find_V1(contigseq, genus_list, V1_dict):
    default_return = "undefined_start"
    for test_genus in genus_list:
        REF_V1 = V1_dict[test_genus]
        if re.search(REF_V1, contigseq):
            return (test_genus)

    return(default_return)

def find_V9(contigseq, genus_list, V9_dict):
    default_return = "undefined_end"
    for test_genus in genus_list:
        REF_V9 = V9_dict[test_genus]
        if re.search(REF_V9, contigseq):
            return (test_genus)
    return(default_return)

def create_pairmatches(original_name, filename):
    # Open the fasta file and an output file
    try:
        input_file = open(filename,'r')
        # make a dictionary using genus as key and the value is the V1 or V9 anchor seq usesd for matching to a contig
        contig_V1_dict = {}
        contig_V1_dict["Bacillus"] = Bacillus_V1
        contig_V1_dict["Enterococcus"] = Enterococcus_V1
        contig_V1_dict["Escherichia"] = Escherichia_V1
        contig_V1_dict["Lactobacillus"] = Lactobacillus_V1
        contig_V1_dict["Listeria"] = Listeria_V1
        contig_V1_dict["Pseudomonas"] = Pseudomonas_V1
        contig_V1_dict["Salmonella"] = Salmonella_V1
        contig_V1_dict["Staphylococcus"] = Staphylococcus_V1

        contig_V9_dict = {}
        contig_V9_dict["Bacillus"] = Bacillus_V9
        contig_V9_dict["Enterococcus"] = Enterococcus_V9
        contig_V9_dict["Escherichia"] = Escherichia_V9
        contig_V9_dict["Lactobacillus"] = Lactobacillus_V9
        contig_V9_dict["Listeria"] = Listeria_V9
        contig_V9_dict["Pseudomonas"] = Pseudomonas_V9
        contig_V9_dict["Salmonella"] = Salmonella_V9
        contig_V9_dict["Staphylococcus"] = Staphylococcus_V9
    
        # list of the keys for matching ends by genus name
        contig_genus_list = []
        contig_genus_list.append("Bacillus")
        contig_genus_list.append("Enterococcus")
        contig_genus_list.append("Escherichia")
        contig_genus_list.append("Lactobacillus")
        contig_genus_list.append("Listeria")
        contig_genus_list.append("Pseudomonas")
        contig_genus_list.append("Salmonella")
        contig_genus_list.append("Staphylococcus")       

        # counting metrics for the entire contig file
        total_contigcount = 0
        trunc_contig_count = 0
        true_count = 0
        chimera_count = 0
        undefined_count = 0
        # for data frame building
        contig_name_col = []
        contig_length_col = []
        left_end_col = []
        right_end_col = []
        # this is a flag value that defines what kind of contig this is: complete (1), chimeric (0), trunc (-1), etc.
        is_complete_col = []
        for each_line in input_file:
            # Is this the sequence line? That occurs in the 2nd line of every 2-line block.
            contig_line = input_file.readline()
            total_contigcount +=1
            # now have the two relevant lines of data.
            left_end_genus = find_V1(contig_line, contig_genus_list, contig_V1_dict)
            right_end_genus = find_V9(contig_line, contig_genus_list, contig_V9_dict)
            # add to the list of contig metrics and outcome for making a data frame col later
            contig_name_col.append(each_line)
            contig_length_col.append(len(contig_line))
            left_end_col.append(left_end_genus)
            right_end_col.append(right_end_genus)
            # only full length contigs can be checked for chimeric, complete, or undefined
            if (len(contig_line) > 1500):
                # if both anchors did not match, this is an undefined contig
                if(left_end_genus == "undefined_start" and right_end_genus == "undefined_end"):
                    undefined_count += 1
                    is_complete_col.append(2)
                else:    
                    if(left_end_genus == right_end_genus):
                        # complete
                        true_count += 1
                        is_complete_col.append(1)
                    if(left_end_genus != right_end_genus):
                        # chimeric
                        chimera_count += 1
                        is_complete_col.append(0)
            else:
                # a truncated contig
                trunc_contig_count +=1
                is_complete_col.append(-1)

        print("total_contigcount: " + str(total_contigcount))
        print("trunc_contig_count: " + str(trunc_contig_count))
        print("true_count: " + str(true_count))
        print("chimera_count: " + str(chimera_count))
        print("undefined_count: " + str(undefined_count))

        input_file.close()

        df = pd.DataFrame({'Name': contig_name_col, 'Length': contig_length_col, 'Start': left_end_col, 'End': right_end_col, 'is_Complete': is_complete_col})
        # save the data frame
        file_name = original_name[:-6] + "_results.csv"
        df.to_csv(file_name)
        # save the count metrics with the original_name file as text
        file_name = original_name[:-6] + "_metrics.txt"
        with open(file_name, 'w') as metrics:
            # Formats the variable into a pretty table string
            metrics.write("total_contigcount: " + str(total_contigcount) + "\n")
            metrics.write("trunc_contig_count: " + str(trunc_contig_count) + "\n")
            metrics.write("true_count: " + str(true_count) + "\n")
            metrics.write("chimera_count: " + str(chimera_count) + "\n")
            metrics.write("undefined_count: " + str(undefined_count) + "\n")
        # now make plots and save them using the original_name 
        bar_plot_maker(original_name, df, contig_genus_list)

    except Exception as e:
        print ("An error was generated.")
        exit ()

def bar_plot_maker(original_name, data_frame, genus_list):
    # this function creates the "nice" looking box plots for counts and percentages for the zymo genus and the discovered chimeras, etc.
    category_list = genus_list.copy()
    category_list.append("chimeric")
    category_list.append("truncated")
    category_list.append("undefined")
    # filter the df based on the type of contig that was discovered by the hunter function
    filtered_complete = data_frame[data_frame['is_Complete'] == 1]
    filtered_chimeric = data_frame[data_frame['is_Complete'] == 0]
    filtered_truncated = data_frame[data_frame['is_Complete'] == -1]
    filtered_undefined = data_frame[data_frame['is_Complete'] == 2]
    total_contigs = len(data_frame)
    # list to store the total counts within each genus and each filtered group that was a chimera, trunc, etc.
    count_list = []
    for genus in genus_list:
        count_list.append((filtered_complete['Start'] == genus).sum())
    # number of "erroneous" contigs in each filter group
    count_list.append(len(filtered_chimeric))
    count_list.append(len(filtered_truncated))
    count_list.append(len(filtered_undefined))
    # color based on normal or erroneous counts
    highlight = {'chimeric', 'truncated','undefined'}
    colors = ['red' if cat in highlight else 'purple' for cat in category_list]
    plt.bar(category_list, count_list, color = colors)
    plt.xlabel('Categories')
    plt.ylabel('Counts')
    plt.title('Complete contigs and chimeras as counts')
    plt.xticks(rotation = 45, ha = 'right')
    plt.tight_layout()
    save_name = original_name[:-6] + '_plot_by_counts.png'
    plt.savefig(save_name, dpi = 300)
    # clear the plot space.
    plt.clf()

    # color based on normal or erroneous counts adjusted as percent of total contigs
    colors = ['orange' if cat in highlight else 'violet' for cat in category_list]
    # convert counts to percentages using list comprehension
    percent_list = [x / total_contigs * 100 for x in count_list]
    plt.bar(category_list, percent_list, color = colors)
    plt.xlabel('Categories')
    plt.ylabel('Percent')
    plt.title('Complete contigs and chimeras as percent')
    plt.xticks(rotation = 45, ha = 'right')
    #plt.tight_layout() already been applied
    save_name = original_name[:-6] + '_plot_by_percent.png'
    plt.savefig(save_name, dpi = 300)
    #plt.show()
    plt.clf()

# -------------------------------------------
#
# MAIN LOOP

if (len(sys.argv) < 2):
        print ("Incorrect number of input files. This program requires 1 input source name.")
        exit ()

chosenfilename = sys.argv[1]
# run FASTA fixwrap first to make it a true 2-line file and feed into the pair matcher
tempfile = FASTA_fixwrap(chosenfilename)
create_pairmatches(chosenfilename, tempfile)
# Remove tempfile
os.remove(tempfile)
