#!/usr/bin/env python3

f = open("pc.wig", "w")

# initialize genome_change variable as a list constituted by 0 with length = genomelength
genome_length = 3079196
genome_change = [0] * genome_length
# open sam file
# sam_file = open(sys.argv[1], 'r')
sam_file = open("lact.sam", 'r')
# read each line of the file

for line in sam_file:
	# if we want to remove whitespace characters
	line = line.rstrip()
	# if line starts with @ we have to skip this line

	if line.startswith('@'):
		continue

	# creation of a list containing all columns of each row. Row is splitted by tab = \t
	fields = line.split("\t")

	# flag indicates that both reads align correctly
	# if ((fields[1] & 3) == 3) and (fields[8] > 0): # => bitwise
	# convertion of flag from integer to binary
	flag = bin(int(fields[1]))

	# get start position and tlen value
	starting_mate_position = int(fields[3])
	mate_length = int(fields[8])

	# if (interesting_flag == '11') and (mate_length > 0):
	if flag.endswith('11') and (mate_length > 0):
		# increment start position by one
		genome_change[starting_mate_position] += 1

		# decrement end position by one
		genome_change[starting_mate_position + mate_length] -= 1

# print genomic profile as a wiggle file
sam_file.close()  # SAM file has been fully read
f.write("fixedStep chrom=genome start=1 step=1 span=1 \n")
current_coverage = 0

# cicle over all positions of the genome
for position in range(genome_length):
	current_coverage += genome_change[position]
	f.write(str(current_coverage) + '\n')

f.close()
