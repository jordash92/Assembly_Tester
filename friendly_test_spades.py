from friendly_fragmenter import *
import os
import time
import Bio

# Assemble fragments in SPAdes

def run_spades(cov, k, f, c):
    os.system("python spades.py --s1 {} -k {} --only-assembler --cov-cutoff {} --phred-offset 64 -o spades_{}".format(f, k, c, cov))

# Counting the number of contigs which appear as an exact match in the genome

def spades_mismatch(contigs, genome):
    count3 = 0
    count4 = 0
    for contig in contigs:
        count4 += 1
        if genome.find(str((contig.seq))) == -1:
            count3 += 1
    return count3, count4

print "Assembling fragments using SPAdes..."
print "May take up to 3 minutes"
print "============================================="

# Call SPAdes and time how long it takes to run in seconds

start_spades_time = time.time()
run_spades(my_coverage, k_value, frag_file, cov_cutoff)
spades_duration = time.time() - start_spades_time

# Assign variables to file containing contigs assembled using SPAdes and to the whole genome 

spades_contigs = Bio.SeqIO.parse("spades_{}/contigs.fasta".format(my_coverage), "fasta")
my_genome = str(get_record("16271976"))

# Call count of mismatches

count3, count4 = spades_mismatch(spades_contigs, my_genome)

print "=============================================="
print "Statistics of SPAdes assembly with", my_coverage, "x coverage"
print "SPAdes took {0:.3f} seconds to asemble the runs".format(spades_duration)
print "A total of {} contigs were assembled".format(count4)
print "{} mismatches out of {} contigs".format(count3, count4)
print "=============================================="

