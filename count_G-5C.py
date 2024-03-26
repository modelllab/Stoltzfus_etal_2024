# Determine fraction of G-5C mutant reads from all fastq files in a directory 
# Note: "fq" files are not actually fqs. They only contain the sequence line from the original fastq files (one per line). 

import glob

#for fq in glob.glob("*.fastq"):
#	print(fq)

# Open output file and write column headers

output = open('counts_G-5C.txt', 'w')
output.write("Sample" + "\t" + "WT_reads" + "\t" + "mut_reads" + "\t" + "other_reads" + "\t" + "total_reads" + "\t" + "fraction_mut" + "\n")

mut = "G"
WT = "C"

print("Time for the loop")

# Iterates through all fastq files in directory you run program from

for fq in glob.glob("*.fastq"):
	print(fq)
	total_reads = 0
	mut_reads = 0
	WT_reads = 0
	other_reads = 0
	short = 0
	
	# Iterates through each line and counts number of mutant/WT reads in -5 position
	file = open(fq)
	for line in file:
		#print(line)
		total_reads += 1
		if len(line) >= 77:
			#print("Long enough")
			base = line[76]
			print(base)
			if base == mut:
					mut_reads += 1
					#print("mutant")
			elif base == WT:
    				WT_reads += 1
    				#print("WT")
			else:
					other_reads += 1
					#print("other")
		else:
			short += 1
			print("Not long enough")
	
	# Calculates mutant fraction and writes results into a new line for each fq file
			
	fraction_mut = float(mut_reads/total_reads)
	output.write(str(fq) + "\t" + str(WT_reads) + "\t" + str(mut_reads) + "\t" + str(other_reads) + "\t" + str(total_reads) + "\t" + str(fraction_mut) + "\n")
	
	print(fq)        	      
	print("Total reads: " + str(total_reads))
	print("Mutant reads: " + str(mut_reads))
	print("WT reads: " + str(WT_reads))
	print("Other reads: " + str(other_reads))
	print("Short reads: " + str(short))
	print("Fraction mut: " + str(fraction_mut))
	print("NEXT")




