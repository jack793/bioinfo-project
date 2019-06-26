#!/usr/bin/env python3

import math
import statistics
import sys

from scripts import pybam
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


def phy_coverage_wig(bam_f):
	print("Physical Coverage started...")

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

	print("Generating .wig file")
	# print genomic profile as a wiggle file
	f.write("fixedStep chrom=genome start=1 step=1 span=1 \n")

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()
	print("done!")


def sequence_coverage_wig(bam_f):
	print("Sequence Coverage started...")

	f = open('../wig-tracks/seq_coverage.wig', 'w')

	# initialize genome_change variable as a list constituted by 0 with length = genomelength
	genome_length = 3079196
	genome_change = [0] * genome_length

	for alignment in bam_f:

		# get start position and tlen value
		starting_mate_position = alignment.sam_pos1  # 4th column
		read_lenght = 100  # 9th column

		if alignment.sam_flag & 4 == 0:  # take all mapped reads
			genome_change[starting_mate_position] += 1

			# decrement end position by one
			genome_change[starting_mate_position + read_lenght] -= 1

	# print genomic profile as a wiggle file
	print("Generating .wig file")
	f.write("fixedStep chrom=genome start=1 step=1 span=1 \n")

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()
	print("done!")


def get_genome_stats(bam_f):
	print("Getting Genome stats...")

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
	avg = total_sum / count_values
	std = statistics.stdev(tlen_list)
	print("Standard Deviation:	" + str(std))
	print("\nTotal Reads:\t{0}\nMax Length:\t{1}\nMin Length:\t{2}\nAverage:\t{3}\nStd.error:\t{4}".format(
		str(count_values), str(max_tlen), str(min_tlen), str(avg), str(std / (math.sqrt(count_values)))))


def distribution_wig(bam_f):
	print("ok")


def avg_inserts_wig(bam_f):
	print("Average genome inserts coverage started...")

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

	print("Generating .wig file")
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
			f.write(str(current_sum/current_coverage)+'\n')
		else:
			f.write('0\n')

	f.close()
	print("done!")


def unique_reads_coverage(bam_f):
	print("Uniquely Mapped Reads Physical Coverage started...")

	f = open('../wig-tracks/unique_reads_coverage.wig', 'w')

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

	print("Generating .wig file")
	# print genomic profile as a wiggle file
	f.write("fixedStep chrom=genome start=1 step=1 span=1 \n")

	current_coverage = 0

	# cicle over all positions of the genome
	for position in range(genome_length):
		current_coverage += genome_change[position]
		f.write(str(current_coverage) + '\n')

	f.close()
	print("done!")


def phy_insert_wig(bam_f):
	cont = 0
	tot = 0
	for read in bam_f:
		if not read.is_unmapped and read.template_length > 0:
			cigarLine = read.cigar
			for cigarType, cigarLength in cigarLine:
				if cigarType == 1:
					cont += 1
					tot += read.template_length
					print("{}:{}-{}".format(read.reference_name, read.reference_start, read.reference_end))


############### DEBUG FUNCTION #################


def debug(bam_f):
	for row in bam_f:
		print(row.bam_mapq)


#################### MAIN ######################


if __name__ == '__main__':

	########################## PART 2 ###########################

	# load sorted Lactobacillus bam file (using pybam library):
	# 	refer to a BAM file sorted by genomic position!
	sorted_bam = pybam.read('../data/lact_sorted.bam')

	# 9) Calculate PHYSICAL COVERAGE, creating related wig file
	# phy_coverage_wig(bam_file)

	# 10) Calculate SEQUENCE COVERAGE, creating related wig file
	# sequence_coverage_wig(bam_file)

	# 11) Calculate insert STATS
	# get_insert_stats(bam_file)

	# 12) Calculate AVERAGE INSERTS LENGTH, creating related wig file
	# avg_inserts_wig(bam_file)

	# 13) Calculate exceeding STD DEVIATION PERCENTAGE , creating related wig file
	# distribution_wig(sorted_bam)

	########################## PART 3 ###########################

	# load Uniquely Mapped reads bam file (filtered from original lact_sorted.bam using samtools)
	#
	# from :
	# https://wabi-wiki.scilifelab.se/display/KB/Filter+uniquely+mapped+reads+from+a+BAM+file#FilteruniquelymappedreadsfromaBAMfile-BWA
	# 	samtools view -h -q 1 -F 4 -F 256 DATA/lact_sorted.bam | grep -v XA:Z | grep -v SA:Z |
	# 		samtools view -b - > DATA/unique.bam
	unique_bam = pybam.read('../data/unique.reads.bam')

	# 14) Calculate UNIQUE READS COVERAGE, creating related wig file
	unique_reads_coverage(unique_bam)

	################# Debug zone ####################

	# debug(bam_file)
