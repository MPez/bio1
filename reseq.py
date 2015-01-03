#Resequencing project

import pysam
import re

pattern = re.compile("/[1|2]$")
insertLength = []


def apriBamFile(bamFileName):
    return pysam.AlignmentFile(bamFileName, "rb")


def getQueryName(query):
    return re.sub(pattern, "", query.query_name)


def comparaFile(file1, file2):
    for read1 in file1:
        print("before for 2")
        print(read1.query_name)
        for read2 in file2:
            print("before if")
            print(read2.query_name)
            if getQueryName(read1) == getQueryName(read2):
                insertLength.append(
                    read1.reference_start - read2.reference_start)
                print(insertLength[len(insertLength)-1])
                break
        print("after for 2")
        print(read1.query_name)
        print(read2.query_name)
    print(read1.query_name)
    print(read2.query_name)

if __name__ == "__main__":
    bamFileName1 = "pass_bam/pass_reads1_sorted_name.bam"
    bamFileName2 = "pass_bam/pass_reads2_sorted_name.bam"

    bamFile1 = apriBamFile(bamFileName1)
    bamFile2 = apriBamFile(bamFileName2)

    comparaFile(bamFile1, bamFile2)

    bamFile1.close()
    bamFile2.close()
