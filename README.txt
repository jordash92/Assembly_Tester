* * * * * Assembly Tester Tutorial * * * * *

Intended purpose of application:

Assembly Tester allows you to determine the quality of the sequence assemblies generate by Velvet and SPAdes of a genome. 
This is done by fragmenting an already solved bacterial genome into 36bp single-ens, reads replicating that of an 
Illumiona sequencer. The number of reads generated is determined by the sequence coverage determined by the user. This 
can be 40, 60, 80, 100, 120 or 140 times coverage.

At the moment, this program only allows the fragmentation and assembly of the bacterial genome Haemophilus influenzae, 
which is 1.8Mb in length. In the future, this program would have ideally been expanded to include the genome of all model 
organisms. This would allow the user to essentialy carry out a test run of an assembly of their desired organism, prior 
to sequenceing a genome to determine the best coverage to sequnece their genome at, as well as choosing an optimal 
assembler for the sequence assembly.

Other future additions would include additional statistics of each assembly, including:
- Max contig length
- Minimum contig length
- Calculation of N50
- Alignment of contigs to original genome
- % genome coverage
- generation of paired-end reads
- option to write stats to file

System requirements:

The following must be installed for application to run:
- Python 2.7.6
- Velvet 1.2.09
- SPAdes 3.6.0

Installation of SPAdes: When in directory which you with so install SPAdes, type the following in to the command line.

wget http://spades.bioinf.spbau.ru/release3.6.0/SPAdes-3.6.0-Linux.tar.gz
tar -xzf SPAdes-3.6.0-Linux.tar.gz
cd SPAdes-3.6.0-Linux/bin/

It is in SPAdes-3.6.0-Linux/bin/ where you must download and extract the contents of Assembly_Tester.zip in order for 
Assembly Tester to run. Assambly Tester is made up of 3 scripts: friendly_fragmenter.py, friendly_test_velvet.py and 
friendly_test_spades.py. friednly_fragmenter.py fragments the genome. friendly_test_velvet.py assembles the fragments 
in velvet. friendly_test_spades.py assembles the fragment using SPAdes and also generates the statistics. Under any 
circumstances, only friendly_test_spades.py should be run, as the functions of the two other scripts are imported into 
this one.

EXAMPLE OF USE / WALKTHROUGH:

EXAMPLE 1: Entry of correct coverage

Step 1. When in SPAdes-3.6.0-Linux/bin/ run friendly_test_spades.py from the command line.

>>> python friendly_test_spades.py

Step 2. You will be requested to enter a genome coverage desired for what would have been the sequenceing step. Coverage 
is the number of reads which a particular base should appear in, given its location within the genome. Options of coverage 
are 40, 60, 80, 100, 120, 140x. For this example, 40x coverage is being selected.

>>> 40

Step 3. Now you sit back and wait for the genome to be fragmented and assembled, firstly in Velvet and then in SPAdes. 
The following should appear on the screen. It describes the parameters which will be used for the assembly in velvet. 
The parameters are specific to the coverage entered. The fasta file in which the fragments are stored is also stated.

Your desired genome coverage is 40x.
The genome of Haemophilus influenzae will now be cut into 2 million fragments.
========================================================== 
Parameters for the assembly:
reads: 2000000
k-value for assembly: 19
coverage cutoff: 18
expected coverage: 19
Haemophilus influenzae has been sucessfully fragmented!
Fragments have been saved to: friendly_fragments_40.fasta
========================================================== 

Step 3. Upon completion of the assembly by both Velvet and SPAdes, the statistics fo the assembly are printed on screen. 
With a coverage of 40x, the results were as follows:

========================================================== 
Statistics of Velvet assembly with 40x coverage
Velvet took 109.435 seconds to asemble the runs
A total of 641 contigs were assembled
108 mismatches out of 641 contigs
==========================================================
Statistics of SPAdes assembly with 40x coverage
SPAdes took 70.688 seconds to asemble the runs
A total of 339 contigs were assembled
181 mismatches out of 339 contigs
==========================================================

NOTE: The results will be different each time due to the generation of random fragments. The statistics contain the duration
of eacvh assembler to run. How many contigs were assembled and finally how many of the contigs assembled did not appear as
an exact match in the original genome.


EXAMPLE 2: Entry of incorrect entry / option to quit

Step 1. When in SPAdes-3.6.0-Linux/bin/ run friendly_test_spades.py from the command line.

>>> python friendly_test_spades.py

Step 2. Upon being asked for a coverage, you enter one which is incorrect (i.e not 40, 60, 80, 100, 120 or 140x), in this 
example 45.

>>> 45

Step 3. You will be asked if you would like to enter a new coverage, which will take you back to the home screen. You would 
do this by pressing 1. Alternatively you can quit Assembly Tester by pressing 2.

>>> 1	# Allows you to enetr new coverage

>>> 2	# Closes Assembly Tester
