"""Insert bacterial genome genbank accession, and extract sequence
data from record. In this case, that of Haemophilus influenzae which
has a genome size of 1.83 Mb."""

from Bio import Entrez
from Bio import SeqIO
import random
import sys

Entrez.email = "j.ashworth2@ncl.ac.uk"

def get_record(id_number):
    """Uses Entrez to parse Genbank online on entry of appropriate Genbank
    id "16271976". Retrieves record of Genbank id. Sequence is extracted
    from record in fasta format and returned."""
    handle = Entrez.efetch(db="nucleotide", rettype="fasta",
                           retmode="text", id="16271976")
    seq_record = SeqIO.read(handle, "fasta")
    handle.close()
    return seq_record.seq

def random_fragment(genome, length):
    """Genome is the .seq attribute of a Genbank record of a whole
    genome. Length is an integer of desired fragment length. A random
    integer between zero and the length of genome minus the size of
    fragment generating to avoid index out of range errors."""
    a = random.randint(0,len(genome)-length)
    # returns slice of genome of sesired length from base "a"
    return genome[a:a+length]

def fragment(file, samples, id):
    """File is a fasta file which fragemtns will be saved to, samples is
    an integer (number of fragments to be generated) and id is a Genbank
    id of genome. For the number of samples stated, random fragments are
    generated of size 36 and written to a fasta file."""
    my_genome = str(get_record(id))
    my_file = open(file, "w")
    # For the number of samples associated woth desired coverage
    for i in range(samples):
        # Give each fragment an individual number, fasta ID and write it to file
        my_file.write(">fragment{}\n".format(i+1))
        # Random fragment from genome of size 36 written to file
        my_file.write(str(random_fragment(my_genome, 36)) + "\n")
    my_file.close()
    return my_genome

def parameters(coverage):
    """Coverage is a integer stated by user. Each coverage has associated
    integers stored in a dictionary: reads, k-value, coverage cutoff
    and expected coverage. Appropriate values returned upon entry of
    coverage. Fasta file is alos created , named according to
    coverage. If user enters unsuitable coverage, given option to
    enter new coverage or quit script.
    """
    reference = {"40": (2000000, 19, 18, 19),
                 "60": (3000000, 29, 20, 21),
                 "80": (4000000, 29, 31, 35),
                 "100": (5000000, 29, 20, 21),
                 "120": (6000000, 31, 18, 19),
                 "140": (7000000, 31, 21, 22)}
    #If incorrect coverage stated by user, asked to enter new coverage or exit
    while coverage not in reference:
        print "ERROR! You have entered an unsuitable coverage, please try again"
        user_retry = raw_input("Would you like to enter a new coverage [1] or leave Assembly Tester[2]?")
        if user_retry == "1":
            # If userselects [1] and enters a new coverage, the old coverage is
            # replaced with new coverage
            coverage = raw_input('Please enter desired genome coverage, this can be 40, 60, 80, 100, 120 or 140 times. For example, if you would like 100x coverage to be generated type 100 and press enter: ')
        # Terminate script if user selects [2] and chooses to exit program
        elif user_retry == "2":
            sys.exit()
        # User did not enter [1] or [2] so asks to try again
        else:
            print "Try again"
    # For coverage selected by user, associate appropriate values in
    # dictonary to following vairables:
    reads, k_value, cov_cutoff, exp_cov = reference[coverage]
    # Name file containing fragments in association to coverage set by user
    frag_file = 'friendly_fragments_{}.fasta'.format(coverage)
    print dots, "\nYour desired genome coverage is {}x.\nThe genome of Haemophilus influenzae will now be cut into {} million fragments.".format(coverage, str(reads)[0])
    return frag_file, reads, k_value, cov_cutoff, exp_cov

dots = "=========================================================="

print dots, "\n * * * * * * * Welcome to Assembly Tester! * * * * * * *\n", dots

"""Ask user for desired genome coverage (40, 60, 80, 100, 120,
140x). Waits for raw input from user before continuing."""

my_coverage = raw_input('Please enter desired genome coverage, this'
'can be 40, 60, 80, 100, 120 or 140 times. For example, if you would'
'like 100x coverage to be generated type 100 and press enter: ')

"""Call parameters to select those appropriate from dictionary for
desired coverage (calculated using Velvet Advisor)."""

frag_file, reads, k_value, cov_cutoff, exp_cov = parameters(my_coverage)

"""Call to fragment the genomeand tho put the files in a "frag_file",
into a certain number of reads and for a particular genbank id, in
this case, that of Haemophilus influenzae."""

fragment(frag_file, reads, "16271976")

print dots, "\nParameters for the assembly:"
print "reads:", reads
print "k-value for assembly:", k_value
print "coverage cutoff:", cov_cutoff
print "expected coverage:", exp_cov
print "Haemophilus influenzae has been sucessfully fragmented!\nFragments have been saved to:", frag_file
