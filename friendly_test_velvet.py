from friendly_fragmenter import *
import os
import time
import Bio

# Assemble fragments in Velvet

def run_velvet(cov, k, f, c, e):
    print "Running Velveth"
    print "=============================================="
    os.system("velveth velvet_{} {} -fasta -short {}".format(cov, k, f))
    print "=============================================="
    print "Running Velvetg"
    print "=============================================="
    os.system("velvetg velvet_{}".format(cov))
    print "=============================================="
    print "Running Velvetg with parameters"
    print "=============================================="
    os.system("velvetg velvet_{} -cov_cutoff {} -exp_cov {} -amos_file yes".format(cov, c, e))

# Counting the number of contigs which appear as an exact match in the genome

def count_mismatch(contigs, genome):
    count = 0
    count2 = 0
    #count_len = []
    for contig in contigs:
        count2 += 1
        if genome.find(str((contig.seq))) == -1:
            count += 1
       # if genome.find(str((contig.seq))) == 1:
    return count, count2

print "Assembling fragments using Velvet..."
print "============================================="

# Call Velvet and time how long it takes to run in seconds

start_velvet_time = time.time()
run_velvet(my_coverage, k_value, frag_file, cov_cutoff, exp_cov)
velvet_duration = time.time() - start_velvet_time

print "Fragments have been assembled!"

# Assign variables to file containing contigs assembled using Velvet and to the whole genome 

contigs = Bio.SeqIO.parse("velvet_{}/contigs.fa".format(my_coverage), "fasta")
genome = str(get_record("16271976"))

# Call count of mismatches

count, count2 = count_mismatch(contigs, genome)

print "=============================================="
print "Statistics of Velvet assembly with", my_coverage, "x coverage"
print "Velvet took {0:.3f} seconds to asemble the runs".format(velvet_duration)
print "A total of {} contigs were assembled".format(count2)
print "{} mismatches out of {} contigs".format(count, count2)
print "=============================================="

