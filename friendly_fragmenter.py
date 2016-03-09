# Insert bacterial genome genbank accession, and extract sequence data from record. In this case, that of Haemophilus influenzae which has a genome size of 1.83 Mb

from Bio import Entrez
from Bio import SeqIO
import random
Entrez.email = "j.ashworth2@ncl.ac.uk"

def get_record(id_number):
    handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id="16271976")
    seq_record = SeqIO.read(handle, "fasta")
    handle.close()
    #print ("%s with %i features and length of %s bp" % (seq_record.id, len(seq_record.features), len(seq_record.seq)))
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

# Fragment genome: put fragments in a file called fragments.fasta, create random 5000000 fragments (will give x100 coverage) of genome associated with genbank id 16271976
# Coverage cutoff and expected coverage were calculated using R, while the k value was calcualted using Velvet Advisor

def parameters(coverage):
    if coverage == "40":
        reads = 2000000
        k_value = 19
        frag_file = 'friendly_fragments_40.fasta'
        cov_cutoff = 18
        exp_cov = 19
        print "Your desired coverage is 40. The genome will now be cut into 2 million fragments."
    elif coverage == "60":
        reads = 3000000
        k_value = 25
        frag_file = 'friendly_fragments_60.fasta'
        cov_cutoff = 18
        exp_cov = 19
        print "Your desired coverage is 60. The genome will now be cut into 3 million fragments."
    elif coverage == "80":
        reads = 4000000
        k_value = 29
        frag_file = 'friendly_fragments_80.fasta'
        cov_cutoff = 31
        exp_cov = 35
        print "Your desired coverage is 80. The genome will now be cut into 4 million fragments."
    elif coverage == "100":
        reads = 5000000
        k_value = 29
        frag_file = 'friendly_fragments_100.fasta'
        cov_cutoff = 20
        exp_cov = 21
        print "Your desired coverage is 100. The genome will now be cut into 5 million fragments."
    elif coverage == "120":
        reads = 6000000
        k_value = 31
        frag_file = 'friendly_fragments_120.fasta'
        cov_cutoff = 18
        exp_cov = 19
        print "Your desired coverage is 120. The genome will now be cut into 6 million fragments."
    elif coverage == "140":
        reads = 7000000
        k_value = 31
        frag_file = 'friendly_fragments_140.fasta'
        cov_cutoff = 21
        exp_cov = 22
        print "Your desired coverage is 140. The genome will now be cut into 7 million fragments."
    else:
        print "ERROR! You have entered an unsuitable coverage, please try again:"
        my_coverage = raw_input('Please enter desired genome coverage, this can be 40, 60, 80, 100, 120 or 140 times. For example, if you would like 100x coverage to be generated type 100 and press enter: ')
        parameters(my_coverage)
    return frag_file, reads, k_value, cov_cutoff, exp_cov

#if __name__ == '__main__':

my_coverage = raw_input('Please enter desired genome coverage, this can be 40, 60, 80, 100, 120 or 140 times. For example, if you would like 100x coverage to be generated type 100 and press enter: ')

frag_file, reads, k_value, cov_cutoff, exp_cov = parameters(my_coverage)

fragment(frag_file, reads, "16271976")

print "============================================="
print "reads: ", reads
print "k-value for assembly: ", k_value
print "coverage cutoff: ", cov_cutoff
print "expected coverage: ", exp_cov
print "Genome has been fragmented! Fragments have been saved to: ", frag_file
