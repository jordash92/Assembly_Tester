from friendly_fragmenter import *
import os
import time
import Bio

'''Assemble fragments in Velvet'''

def run_velvet(cov, k, f, c, e):
    """Runs Velvet from command line with appropriate parameters of which
    there are three steps: velveth, velvetg and velvetg with parameters.
    Where cov is coverage, k is k-value, f is file containing fragments,
    c is coverage cutoff and e is expected coverage. "-amos_file yes"
    creates a file allowing you to view assembled fragments in Tablet.
    Results of assembly are stored in file names after the assembler and
    coverage."""
    print "Running Velveth\n", dots
    os.system("velveth velvet_{} {} -fasta -short {}".format(cov, k, f))
    print dots, "\nRunning Velvetg\n", dots
    os.system("velvetg velvet_{}".format(cov))
    print dots, "\nRunning Velvetg with parameters\n", dots
    os.system("velvetg velvet_{} -cov_cutoff {} -exp_cov {} -amos_file yes"
              .format(cov, c, e))

def count_mismatch(contigs, genome):
    """Count number of contigs assembled by Velvet and how many of those
    in original the genome. Contigs is a fasta file, while the genome is """
    count = 0
    count2 = 0
    # count the number of contigs which Velvet has assembled
    for contig in contigs:
        count2 += 1
        """Count the number of contigs which Velvet has assembled which DO NOT
        appear as an exact match in the original genome"""
        if genome.find(str((contig.seq))) == -1:
            count += 1
    return count, count2

# Call Velvet and time how long it takes for the assembly to run in seconds

print dots, "\nAssembling fragments using Velvet...\n", dots
start_velvet_time = time.time()
run_velvet(my_coverage, k_value, frag_file, cov_cutoff, exp_cov)
velvet_duration = time.time() - start_velvet_time
print dots, "\nVelvet has sucessfully assembled the fragments!\n", dots

"""Assign variables to file containing contigs assembled using Velvet
and to the whole genome."""

velvet_contigs = Bio.SeqIO.parse("velvet_{}/contigs.fa"
                                 .format(my_coverage), "fasta")
genome = str(get_record("16271976"))

# Call count of mismatches generated by Velvet

count, count2 = count_mismatch(velvet_contigs, genome)
