# Determine fraction of reads with WT (GGG) and mutant (GGA) PAM sequences from competition experiments. 
# After demultiplexing reads, run this program from the directory containing resulting "fq" files
	# Note: "fq" files are not fastqs. They are text files containing just the sequence lines from the original fastqs. 

# Import glob, allowing us to iterate through all fq files within directory.
import glob

# Create output file and write column headers.
output = open('counts_PAM.txt', 'w')
output.write("Sample" + "\t" + "WT_reads" + "\t" + "mut_reads" + "\t" + "other_reads" + "\t" + "total_reads" + "\t" + "fraction_mut" + "\n")

# Define WT and mutant base identities (Note: PAM is the on reverse strand, so the WT base will be C, and the mutant base will be T).
mut = "T"
WT = "C"

# Iterate through all fq files in directory from which you run the program.
for fq in glob.glob("*.fq"):
	print(fq)
	total_reads = 0
	mut_reads = 0
	WT_reads = 0
	other_reads = 0
	short = 0
	# For each fq file, iterate through each line and count number of mutant/WT reads in -5 position.
	file = open(fq)
	for line in file:
		total_reads += 1
		# If a read is long enough, determine if the sequence is WT or mutant.
		if len(line) >= 69:
			PAM = line[69]
			if PAM == mut:
					mut_reads += 1
			elif PAM == WT:
    				WT_reads += 1
			else:
					other_reads += 1
		# If a read is not long enough, do not count it.
		else:
			short += 1
	
	# Calculate mutant fraction and write results into a new line for each fq file.
	fraction_mut = float(mut_reads/total_reads)
	output.write(str(fq) + "\t" + str(WT_reads) + "\t" + str(mut_reads) + "\t" + str(other_reads) + "\t" + str(total_reads) + "\t" + str(fraction_mut) + "\n")
	




