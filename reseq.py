#Resequencing project

import pysam
import re

def apriBamFile(bamFileName):
	return pysam.AlignmentFile(bamFileName, "rb")

def comparaFile(file1, file2):
	for read1 in bamFile1:
		for read2 in bamFile2:
			if (re.sub("/[1|2]$", "", read1.query_name) == re.sub("/[1|2]$", "", read2.query_name))
				print ("equal")
			else:
				print("not equal")

if __name__ == "__main__":
	bamFileName1 = "pass_bam/pass_reads1_sorted_name.bam"
	bamFileName2 = "pass_bam/pass_reads2_sorted_name.bam"

	bamFile1 = apriBamFile(bamFileName1)
	bamFile2 = apriBamFile(bamFileName2)

	comparaFile(bamFile1,bamFile2)