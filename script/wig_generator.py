from script import pybam

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


for alignment in pybam.read('data/lact_sorted.bam'):
        # get flag,
        print(alignment.sam_flag, alignment.sam_pos1, alignment.sam_tlen)


# tupla: (quante pairs son presenti, lunghezza media della somma dei pairs )


