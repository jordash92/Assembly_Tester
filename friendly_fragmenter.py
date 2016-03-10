# Insert bacterial genome genbank accession, and extract sequence data from record. In this case, that of Haemophilus influenzae which has a genome size of 1.83 Mb

from Bio import Entrez
from Bio import SeqIO
import random
import sys

Entrez.email = "j.ashworth2@ncl.ac.uk"



def get_record(id_number):
    handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="16271976")
    seq_record = SeqIO.read(handle, "fasta")
    handle.close()
    return seq_record.seq

def random_fragment(genome, length):
    a = random.randint(0,len(genome)-length)
    return genome[a:a+length]

def fragment(file, samples, id):
    my_genome = str(get_record(id))
    my_file = open(file, "w")
    for i in range(samples):
        my_file.write(">fragment{}\n".format(i+1))
        my_file.write(str(random_fragment(my_genome, 36)) + "\n")
    my_file.close()
    return my_genome

# Parameters for the assembly chosen depending on coverage stated by user
# Coverage cutoff and expected coverage were calculated using R, while the k value was calcualted using Velvet Advisor

def parameters(coverage):
    reference = {"40": (2000000, 19, 18, 19),
                 "60": (3000000, 29, 20, 21),
                 "80": (4000000, 29, 31, 35),
                 "100": (5000000, 29, 20, 21),
                 "120": (6000000, 31, 18, 19),
                 "140": (7000000, 31, 21, 22)}
    # if incorrect coverage stated by user, asked to enter new coverage or exit program
    while coverage not in reference:
        print "ERROR! You have entered an unsuitable coverage, please try again"
        user_retry = raw_input("Would you like to enter a new coverage [1] or leave Assembly Tester [2]?")
        if user_retry == "1":
            # if user chooses to enter new coverage, replaces old coverage for new coverage
            coverage = raw_input('Please enter desired genome coverage, this can be 40, 60, 80, 100, 120 or 140 times. For example, if you would like 100x coverage to be generated type 100 and press enter: ')
        # terminate script if user chooses to exit program
        elif user_retry == "2":
            sys.exit()
        # user did not enter 1 or 2 so asks to try again
        else:
            print "Try again"
    # for coverage selected by user, associate appropriate values in dictonary to following vairables:
    reads, k_value, cov_cutoff, exp_cov = reference[coverage]
    # name file containing fragments is association to coverage set by user
    frag_file = 'friendly_fragments_{}.fasta'.format(coverage)
    print dots, "\nYour desired coverage is {}x.\nThe genome of Haemophilus influenzae will now be cut into {} million fragments.".format(coverage, str(reads)[0])
    return frag_file, reads, k_value, cov_cutoff, exp_cov

dots = "=========================================================="

print dots, "\n * * * * * * * Welcome to Assembly Tester! * * * * * * *\n", dots

# Ask user for desired coverage. Waits for raw input from user before continuing.

my_coverage = raw_input('Please enter desired genome coverage, this can be 40, 60, 80, 100, 120 or 140 times. For example, if you would like 100x coverage to be generated type 100 and press enter: ')

# Call parameters to select those appropriate from dictionary for desired coverage (calculated using Velvet Advisor)

frag_file, reads, k_value, cov_cutoff, exp_cov = parameters(my_coverage)

# Call fragment genome

fragment(frag_file, reads, "16271976")

print dots, "\nParameters for the assembly:"
print "reads:", reads
print "k-value for assembly:", k_value
print "coverage cutoff:", cov_cutoff
print "expected coverage:", exp_cov
print "Haemophilus influenzae has been sucessfully fragmented!\nFragments have been saved to:", frag_file


