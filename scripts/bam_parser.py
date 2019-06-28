#!/usr/bin/env python3

import math
import statistics
import sys
import os

import pybam


# https://github.com/luidale/pybam [Currently only Dynamic Parsing works]


# bam --------------------- All the bytes that make up the current alignment ("read"),
#                                   still in binary just as it was in the BAM file. Useful
#                                   when creating a new BAM file of filtered alignments.
#
# ===============================================================================================
#
#       sam --------------------- The current alignment in SAM format.
#       sam_bin ----------------- The bin value of the alignment (used for indexing reads).
#                                   Please refer to section 5.3 of the SAM spec for how this
#                                   value is calculated.
#       sam_block_size ---------- The number of bytes the current alignment takes up in the BAM
#                                   file minus the four bytes used to store the block_size value
#                                   itself. Essentially sam_block_size +4 == bytes needed to store
#                                   the current alignment.
#       sam_cigar_list ---------- A list of tuples with 2 values per tuple:
#                                   the number of bases, and the CIGAR operation applied to those
#                                   bases. Faster to calculate than sam_cigar_string.
#       sam_cigar_string -------- [6th column in SAM] The CIGAR string, as per the SAM format.
#                                   Allowed values are "MIDNSHP=X".
#       sam_flag ---------------- [2nd column in SAM] The FLAG number of the alignment.
#       sam_l_read_name --------- The length of the QNAME plus 1 because the QNAME is terminated
#                                   with a NUL byte.
#       sam_l_seq --------------- The number of bases in the seq. Useful if you just want to know
#                                   how many bases are in the SEQ but do not need to know what those
#                                   bases are (which requires more decoding effort).
#       sam_mapq ---------------- [5th column in SAM] The Mapping Quality of the current alignment.
#       sam_n_cigar_op ---------- The number of CIGAR operations in the CIGAR field. Useful if one
#                                   wants to know how many CIGAR operations there are, but does not
#                                   need to know what they are.
#       sam_next_refID ---------- The sam_refID of the alignment's mate (if any). Note that as per
#                                   sam_refID, this value can be negative and is not the actual
#                                   chromosome name (see sam_pnext1).
#       sam_pnext0 -------------- The 0-based position of the alignment's mate. Note that in SAM all
#                                   positions are 1-based, but in BAM they are stored as 0-based.
#                                   Unlike sam_pnext1, negative values are kept as negative values
#                                   here, essentially giving you the value as it was stored in BAM.
#       sam_pnext1 -------------- [8th column in SAM] The 1-based position of the alignment's mate.
#                                   Note that in SAM format values less than 1 are converted to "0"
#                                   for "no data", and sam_pnext1 will also do this.
#       sam_pos0 ---------------- The 0-based position of the alignment. Note that in SAM all
#                                   positions are 1-based, but in BAM they are stored as 0-based.
#                                   Unlike sam_pos1, negative values are kept as negative values,
#                                   essentially giving one the decoded value as it was stored.
#       sam_pos1 ---------------- [4th column in SAM] The 1-based position of the alignment. Note
#                                   that in SAM format values less than 1 are converted to "0" for
#                                   "no data" and sam_pos1 will also do this.
#       sam_qname --------------- [1st column in SAM] The QNAME (fragment ID) of the alignment.
#       sam_qual ---------------- [11th column in SAM] The QUAL value (quality scores per DNA base
#                                   in SEQ) of the alignment.
#       sam_refID --------------- The chromosome ID (not the same as the name!).
#                                   Chromosome names are stored in the BAM header (file_chromosomes),
#                                   so to convert refIDs to chromsome names one needs to do:
#                                   "my_bam.file_chromosomes[read.sam_refID]" (or use sam_rname)
#                                   But for comparisons, using the refID is much faster that using
#                                   the actual chromosome name (for example, when reading through a
#                                   sorted BAM file and looking for where last_refID != this_refID)
#                                   Note that when negative the alignment is not aligned, and thus one
#                                   must not perform my_bam.file_chromosomes[read.sam_refID]
#                                   without checking that the value is positive first.
#       sam_rname --------------- [3rd column in SAM] The actual chromosome/contig name for the
#                                   alignment. Will return "*" if refID is negative.
#       sam_rnext --------------- [7th column in SAM] The chromosome name of the alignment's mate.
#                                   Value is "*" if unmapped. Note that in a SAM file this value
#                                   is "=" if it is the same as the sam_rname, however pybam will
#                                   only do this if the user prints the whole SAM entry with "sam".
#       sam_seq ----------------- [10th column in SAM] The SEQ value (DNA sequence of the alignment).
#                                   Allowed values are "ACGTMRSVWYHKDBN and =".
#       sam_tags_list ----------- A list of tuples with 3 values per tuple: a two-letter TAG ID, the
#                                   type code used to describe the data in the TAG value (see SAM spec.
#                                   for details), and the value of the TAG. Note that the BAM format
#                                   has type codes like "c" for a number in the range -127 to +127,
#                                   and "C" for a number in the range of 0 to 255.
#                                   In a SAM file however, all numerical codes appear to just be stored
#                                   using "i", which is a number in the range -2147483647 to +2147483647.
#                                   sam_tags_list will therefore return the code used in the BAM file,
#                                   and not "i" for all numbers.
#       sam_tags_string --------- [12th column a SAM] Returns the TAGs in the same format as would be found
#                                   in a SAM file (with all numbers having a signed 32bit code of "i").
#       sam_tlen ---------------- [9th column in SAM] The TLEN value.


################ FUNCTIONS #################

def total_phy_coverage(bam_f):
	print("Physical Coverage started...\n")

	f = open('../wig-tracks/tot_phy_coverage.wig', 'w')

	# initialize genome_change variable as a list constituted by 0 with length = genomelength
	genome_length = 3079196
	genome_change = [0] * genome_length

	for alignment in bam_f:
		# conversion of flag from integer to binary
		flag = bin(alignment.sam_flag)

		# get start position and tlen value
		starting_mate_position = alignment.sam_pos1  # 4th column
		template_length = alignment.sam_tlen  # 9th column

		if flag.endswith('11'):
			# increment start mate position by one
			genome_change[starting_mate_position] += 1
			# decrement end mate position by one
			genome_change[starting_mate_position + template_length] -= 1

	print("Generating .wig file\n")
	# print genomic profile as a wiggle file
	f.write("fixedStep chrom=genome start=1 step=1 span=1 \n")

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()
	print("done!")


def phy_coverage(bam_f):
	print("Physical Coverage started...\n")

	f = open('../wig-tracks/phy_coverage.wig', 'w')

	# initialize genome_change variable as a list constituted by 0 with length = genomelength
	genome_length = 3079196
	genome_change = [0] * genome_length

	for alignment in bam_f:
		# conversion of flag from integer to binary
		flag = bin(alignment.sam_flag)

		# get start position and tlen value
		starting_mate_position = alignment.sam_pos1  # 4th column
		template_length = alignment.sam_tlen  # 9th column

		if 0 < template_length <= 3000:

			# if (interesting_flag == '11') and (mate_length > 0)
			# 	11 -> read paired and read paired in proper pair
			if flag.endswith('11'):
				# increment start position by one
				genome_change[starting_mate_position] += 1

				# decrement end position by one
				genome_change[starting_mate_position + template_length] -= 1

	print("Generating .wig file\n")
	# print genomic profile as a wiggle file
	f.write("fixedStep chrom=genome start=1 step=1 span=1 \n")

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()
	print("done!")

	input("\npress any key to return in the menu...\n")
	return main_menu()  # <----- Recall Main tool view


def sequence_coverage(bam_f):
	print("Sequence Coverage started...\n")

	f = open('../wig-tracks/seq_coverage.wig', 'w')

	# initialize genome_change variable as a list constituted by 0 with length = genomelength
	genome_length = 3079196
	genome_change = [0] * genome_length

	for alignment in bam_f:

		# get start position and tlen value
		starting_mate_position = alignment.sam_pos1  # 4th column
		read_lenght = 100  # 9th column

		if alignment.sam_flag & 4 == 0:  # take all mapped reads: [3rd bit] in bam format == 0
			genome_change[starting_mate_position] += 1

			# decrement end position by one
			genome_change[starting_mate_position + read_lenght] -= 1

	# print genomic profile as a wiggle file
	print("Generating .wig file\n")
	f.write("fixedStep chrom=genome start=1 step=1 span=1 \n")

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()
	print("done!")

	input("\npress any key to return in the menu...\n")
	return main_menu()  # <----- Recall Main tool view


def get_genome_stats(bam_f):
	print("Getting Genome stats...\n")

	# init vars
	count_values = 0
	total_sum = 0
	max_tlen = 0
	min_tlen = sys.maxsize
	tlen_list = []

	for alignment in bam_f:

		# get tlen value
		template_length = alignment.sam_tlen  # 9th column

		if 0 < template_length <= 3000:  # manual filter discarding OUTLIERS

			count_values += 1
			if template_length > max_tlen:
				max_tlen = template_length

			if template_length < min_tlen:
				min_tlen = template_length

			total_sum += template_length
			tlen_list.append(template_length)

	print("done..Calculating statistics..\n\n")
	# Standard deviation calc
	average = total_sum / count_values
	stdev = statistics.stdev(tlen_list)
	print("Standard Deviation:	" + str(stdev))
	print("\nTotal Reads:\t{0}\nMax Length:\t{1}\nMin Length:\t{2}\nAverage:\t{3}\nStd.error:\t{4}".format(
		str(count_values), str(max_tlen), str(min_tlen), str(average), str(stdev / (math.sqrt(count_values)))))

	input("\npress any key to return in the menu...\n")
	return main_menu()  # <----- Recall Main tool view


def two_distribution_percentage(bam_f, average, stdev):
	print("Standard Deviaton Distribution started...\n")

	f = open('../wig-tracks/std_distribution(temp).wig', 'w')

	# initialize genome_change variable as a list constituted by 0 with length = genomelength
	genome_length = 3079196
	genome_change = [0] * genome_length

	for alignment in bam_f:
		# conversion of flag from integer to binary
		flag = bin(alignment.sam_flag)

		# get start position and tlen value
		starting_mate_position = alignment.sam_pos1  # 4th column
		template_length = alignment.sam_tlen  # 9th column

		if flag.endswith('11') and (template_length > 0) and \
			((template_length < average - stdev * 2) or (template_length > average + stdev * 2)):

			# if (interesting_flag == '11') and (mate_length > 0)
			# 	11 -> read paired and read paired in proper pair
			if flag.endswith('11'):
				# increment start position by one
				genome_change[starting_mate_position] += 1

				# decrement end position by one
				genome_change[starting_mate_position + template_length] -= 1

	print("Generating temp .wig file\n")
	# print genomic profile as a wiggle file

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()

	print("calculating percentage coverage distribution in 'std_%distribution.wig'...\n")
	get_percent_match("phy_coverage.wig", "std_distribution(temp).wig", "std_%distribution.wig")

	print("Removing 'std_distribution(temp)' wig file..\n")
	os.remove("../wig-tracks/std_distribution(temp).wig")
	print("	std_distribution(temp).wig removed.\n")

	print("done!")

	input("\npress any key to return in the menu...\n")
	return main_menu()  # <----- Recall Main tool view


def avg_inserts_coverage(bam_f):
	print("Average genome inserts coverage started...\n")

	f = open('../wig-tracks/avg_inserts_coverage.wig', 'w')

	# initialize genome_change variable as a list constituted by 0 with length = genomelength
	genome_length = 3079196
	genome_change = [0] * genome_length
	genome_reads = [0] * genome_length

	for alignment in bam_f:
		# conversion of flag from integer to binary
		flag = bin(alignment.sam_flag)

		# get start position and tlen value
		starting_mate_position = alignment.sam_pos1  # 4th column
		template_length = alignment.sam_tlen  # 9th column

		if 0 < template_length <= 3000:

			# if (interesting_flag == '11') and (mate_length > 0)
			# 	11 -> read paired and read paired in proper pair
			if flag.endswith('11'):
				# increment start position by one
				genome_change[starting_mate_position] += 1
				genome_reads[starting_mate_position] += template_length

				# decrement end position by one
				genome_change[starting_mate_position + template_length] -= 1
				genome_reads[starting_mate_position + template_length] -= template_length

	print("Generating .wig file\n")
	# print genomic profile as a wiggle file
	f.write("fixedStep chrom=genome start=1 step=1 span=1 \n")

	# init useful vars for write wig avg track
	current_coverage = 0
	current_sum = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		current_sum += genome_reads[position]

		if current_coverage > 0:
			f.write(str(current_sum / current_coverage) + '\n')
		else:
			f.write('0\n')

	f.close()
	print("done!")

	input("\npress any key to return in the menu...\n")
	return main_menu()  # <----- Recall Main tool view


def unique_reads_coverage(bam_f):
	print("Uniquely Mapped Reads Physical Coverage started...\n")

	f = open('../wig-tracks/unique_coverage.wig', 'w')

	# initialize genome_change variable as a list constituted by 0 with length = genomelength
	genome_length = 3079196
	genome_change = [0] * genome_length

	for alignment in bam_f:
		# conversion of flag from integer to binary
		flag = bin(alignment.sam_flag)

		# get start position and tlen value
		starting_mate_position = alignment.sam_pos1  # 4th column
		template_length = alignment.sam_tlen  # 9th column

		if 0 < template_length <= 3000:

			# if (interesting_flag == '11') and (mate_length > 0)
			# 	11 -> read paired and read paired in proper pair
			if flag.endswith('11'):
				# increment start position by one
				genome_change[starting_mate_position] += 1

				# decrement end position by one
				genome_change[starting_mate_position + template_length] -= 1

	print("Generating 'unique_coverage.wig' file\n")
	# print genomic profile as a wiggle file

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()

	print("done!")

	input("\npress any key to return in the menu...\n")
	return main_menu()  # <----- Recall Main tool view


def multiple_reads_coverage(bam_f):
	print("Multiple Mapped Reads Coverage started...\n")

	f = open('../wig-tracks/multiple_reads_coverage.wig', 'w')

	# initialize genome_change variable as a list constituted by 0 with length = genomelength
	genome_length = 3079196
	genome_change = [0] * genome_length

	for alignment in bam_f:
		# conversion of flag from integer to binary
		flag = bin(alignment.sam_flag)
		multi_flag = flag[-9:]

		# get start position and tlen value
		starting_mate_position = alignment.sam_pos1  # 4th column
		template_length = alignment.sam_tlen  # 9th column

		if template_length <= 3000 and multi_flag.startswith('1'):

			# if (interesting_flag == '11') and (mate_length > 0)
			# 	11 -> read paired and read paired in proper pair
			if flag.endswith('11'):
				# increment start position by one
				genome_change[starting_mate_position] += 1

				# decrement end position by one
				genome_change[starting_mate_position + template_length] -= 1

	print("Generating 'multiple_reads_coverage.wig' file\n")
	# print genomic profile as a wiggle file

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()
	print("done!")

	input("\npress any key to return in the menu...\n")
	return main_menu()  # <----- Recall Main tool view


def oriented_mates_percentage(bam_f):
	print("Oriented mates percentage calculation started...\n")

	f1 = open('../wig-tracks/out_oriented(temp).wig', 'w')
	f2 = open('../wig-tracks/in_oriented(temp).wig', 'w')
	f3 = open('../wig-tracks/wrong_oriented(temp).wig', 'w')

	# initialize 3 variables as list constituted by 0 with length = genome_length
	genome_length = 3079196

	genome_out = [0] * genome_length
	genome_in = [0] * genome_length
	genome_wrong = [0] * genome_length

	for alignment in bam_f:
		# conversion of flag from integer to binary
		flag = bin(alignment.sam_flag)[-6:]  # INTERESTED flag: take last 6 bits

		# get start position and tlen value
		start_pos = alignment.sam_pos1  # 4th column
		tlen = alignment.sam_tlen  # 9th column

		if tlen <= 3000 and flag.endswith('1'):  # filter + paired

			## case: <-- <-- & --> --> as reads are WRONG MAPPED,
			if flag.startswith('11') or flag.startswith('00'):
				if tlen > 0:
					genome_wrong[start_pos] += 1
					genome_wrong[start_pos + tlen] -= 1
				else:
					genome_wrong[start_pos + tlen + 1] += 1
					genome_wrong[start_pos + 1] -= 1


			## case: <-- --> positive strand and positive length(10 & l>0),
			elif flag.startswith('10') and tlen > 0:
				genome_out[start_pos] += 1
				genome_out[start_pos + tlen] -= 1

			## case: <-- --> negative strand and negative length(01 & l<0),
			elif flag.startswith('01') and tlen < 0:
				genome_out[start_pos + tlen + 1] -= 1
				genome_out[start_pos + 1] += 1


			## case: --> <-- positive strand but negative length(10 & l<0)
			elif flag.startswith('01') and tlen > 0:
				genome_in[start_pos] += 1
				genome_in[start_pos + tlen] -= 1

			## case: --> <-- negative strand but positive length(01 & l>0)
			elif flag.startswith('10') and tlen < 0:
				genome_in[start_pos + tlen + 1] += 1
				genome_in[start_pos + 1] -= 1

	print("Generating 3 temp .wig file...please wait a moment\n")
	# print genomic profile as a wiggle file

	current_coverage1 = 0  # out
	current_coverage2 = 0  # right
	current_coverage3 = 0  # left

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage1 += genome_out[position]
		current_coverage2 += genome_in[position]
		current_coverage3 += genome_wrong[position]

		f1.write(str(current_coverage1) + '\n')
		f2.write(str(current_coverage2) + '\n')
		f3.write(str(current_coverage3) + '\n')

	f1.close()
	f2.close()
	f3.close()
	print("Calculation the percentage match for the three files..please wait\n")

	get_percent_match("phy_coverage.wig", "out_oriented(temp).wig", "out_%oriented.wig")
	get_percent_match("phy_coverage.wig", "in_oriented(temp).wig", "in_%oriented.wig")
	get_percent_match("phy_coverage.wig", "wrong_oriented(temp).wig", "wrong_%oriented.wig")

	print("Removing temp .wig files...\n")

	os.remove("../wig-tracks/out_oriented(temp).wig")
	print("out_oriented(temp).wig removed.")
	os.remove("../wig-tracks/in_oriented(temp).wig")
	print("in_oriented(temp).wig removed.")
	os.remove("../wig-tracks/wrong_oriented(temp).wig")
	print("wrong_oriented(temp).wig removed.\n")

	print("done!")

	input("\npress any key to return in the menu...\n")
	return main_menu()  # <----- Recall Main tool view


def single_mates_percentage(bam_f):
	print("Single Mates percentage calculation started...\n")

	f = open('../wig-tracks/single_mates(temp).wig', 'w')

	# initialize genome_change variable as a list constituted by 0 with length = genomelength
	genome_length = 3079196
	genome_change = [0] * genome_length

	for alignment in bam_f:
		# conversion of flag from integer to binary
		single_flag = bin(alignment.sam_flag)[-4:]  # 4th bit: 0x4 stand for SEGMENT UNMAPPED

		# get start position and tlen value
		start_pos = alignment.sam_pos1  # 4th column: start mate pos

		if single_flag.startswith('01'):  # not mapped reads(1) but complementary is(0)
			# Only start position increase, due to SAM specification that say:
			# 	TLEN value is set as 0 for single-segment template...
			genome_change[start_pos] += 1

	print("Generating .wig file\n")
	# print genomic profile as a wiggle file

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()

	print("Calculate total physical coverage and put it in a wig file...wait a moment..\n")
	total_phy_coverage(bam_f)

	print("Generating 'single_mates_%coverage.wig' file")
	get_percent_match("tot_phy_coverage.wig", "single_mates(temp).wig", "single_mates_%coverage.wig")

	print("Removing temp 'single_mates(temp)' wig files...\n")

	os.remove("../wig-tracks/single_mates(temp).wig")
	print("single_mates(temp).wig removed.\n")

	print("done!")

	input("\npress any key to return in the menu...\n")
	return main_menu()  # <----- Recall Main tool view


############### DEBUG FUNCTIONS #################


def get_percent_match(inp: str, temp: str, out: str):
	with open('../wig-tracks/' + inp, "r") as f1, open('../wig-tracks/' + temp, "r") as f2, \
		open('../wig-tracks/' + out, "w") as out:

		v1 = []
		v2 = []
		i = 0
		wig_String = False

		for line in f1.readlines():
			if not wig_String:
				out.write("fixedStep chrom=genome start=1 step=1 span=1\n")
				wig_String = True
			else:
				line = line.rstrip()
				v1.insert(i, int(line))
				i = i + 1

		i = 0

		for line in f2.readlines():
			line = line.rstrip()
			v2.insert(i, int(line))
			i = i + 1

		for x in range(len(v1)):
			if v1[x] is 0:
				out.write("0.0 \n")
			else:
				out.write(str(100 * v2[x] / v1[x]) + "\n")


def debug(bam_f):
	print("debugging...\n")
	cont = 0
	for row in bam_f:
		lenght = row.sam_tlen
		flag = bin(int(row.sam_flag))
		print(flag, lenght, cont)


#################### MAIN ######################


if __name__ == '__main__':
	########################## PART 2 ###########################

	# load sorted Lactobacillus bam file (using pybam library):
	# 	refer to a BAM file sorted by genomic position!
	sorted_bam = pybam.read('../data/lact_sorted.bam')

	# 9) Calculate PHYSICAL COVERAGE, creating related wig file
	# phy_coverage(sorted_bam)

	# 10) Calculate SEQUENCE COVERAGE, creating related wig file
	# sequence_coverage(sorted_bam)

	# 11) Calculate INSERT STATS
	# get_genome_stats(sorted_bam)

	# 12) Calculate AVERAGE INSERTS LENGTH, creating related wig file
	# avg_inserts_coverage(sorted_bam)

	# Saved values from get_genome_stats() function above
	avg = 2101.0225496051385
	std = 201.66577043621606


	# 13) Calculate exceeding STD DEVIATION PERCENTAGE, creating related wig file
	# distribution_percentage(sorted_bam, avg, std)

	########################## PART 3 ###########################

	# load Uniquely Mapped reads bam file (filtered from original lact_sorted.bam using samtools)
	#
	# from :
	# https://wabi-wiki.scilifelab.se/display/KB/Filter+uniquely+mapped+reads+from+a+BAM+file#FilteruniquelymappedreadsfromaBAMfile-BWA
	# 	samtools view -h -q 1 -F 4 -F 256 DATA/lact_sorted.bam | grep -v XA:Z | grep -v SA:Z |
	# 		samtools view -b - > DATA/unique.bam
	# unique_bam = pybam.read('../data/unique.reads.bam')

	# 14) Calculate UNIQUE READS COVERAGE, creating related wig file
	# unique_reads_coverage(unique_bam)

	# 15) Calculate MULTIPLE READS COVERAGE, creating related wig file
	# multiple_reads_coverage(sorted_bam)

	# 16) Calculate ORIENTED MATES PERCENTAGE, creating related wig file
	# oriented_mates_percentage(sorted_bam)

	# 17) Calculate SINGLE MATES PERCENTEAGE, creating related wig file
	# single_mates_percentage(sorted_bam)

	################# Debug zone ####################

	# debug(sorted_bam)

	#################################################
	#################  PARSER TOOL  #################
	#################################################

	# =======================
	#     MENUS FUNCTIONS
	# =======================

	# Main menu
	def main_menu():
		os.system('clear')

		print("Welcome,\n")
		print("Please choose the functionality you want to start:\n")

		print("######## PART II #########")
		print("[1] Physical coverage")
		print("[2] Sequence coverage")
		print("[3] Genome Inserts stats")
		print("[4] Average inserts length coverage")
		print("[5] Inserts length percentage(%) exceeding ST Deviation\n")

		print("####### PART III ########")
		print("[6] Unique reads coverage")
		print("[7] Multiple reads coverage")
		print("[8] Oriented mates percentage(%)")
		print("[9] Single mates percentage(%)\n")

		print("[0]. Quit")

		choice = input(">>  ")  # Ask for a string in std input
		run_menu(choice)  # Execute user choice


	# Execute menu
	def run_menu(choice):
		os.system('clear')  # Clear OS Terminal CLI
		ch = choice.lower()
		if ch == '':
			print("Please, select a valid option")
			choice = input(">> ")
			run_menu(choice)
		else:
			if int(ch) == 5:
				two_distribution_percentage(sorted_bam, avg, std)  # two_distr.. need 3 args
			elif int(ch) == 0:
				exit_program()
			elif 0 < int(ch) < 10:  # VALID OPTION
				menu_actions[int(ch)](sorted_bam)
			else:  # NOT VALID OPTION
				print("Invalid selection, please try again.\n")
				main_menu()
		return


	# Exit program
	def exit_program():
		print("Exit from program..see you soon")
		sys.exit()


	# =======================
	#    MENUS DEFINITIONS
	# =======================

	menu_actions = {
		1: phy_coverage,
		2: sequence_coverage,
		3: get_genome_stats,
		4: avg_inserts_coverage,
		5: two_distribution_percentage,
		6: unique_reads_coverage,
		7: multiple_reads_coverage,
		8: oriented_mates_percentage,
		9: single_mates_percentage,
		0: exit_program
	}

	# =======================
	#      MAIN PROGRAM
	# =======================

	main_menu()
