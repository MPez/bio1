#Resequencing project

import pysam
import re
import math

bam_file_name1 = "pass_bam/pass_reads1_sorted_name.bam"
bam_file_name2 = "pass_bam/pass_reads2_sorted_name.bam"
insert_length = []


def apri_bam_file(num):
    if num == 1:
        return pysam.AlignmentFile(bam_file_name1, "rb")
    else:
        return pysam.AlignmentFile(bam_file_name2, "rb")


def getQueryName(query):
    qname_pattern = re.compile("/[1|2]$")
    return re.sub(qname_pattern, "", query.query_name)


def skip_read(bam_file, num):
    it = bam_file.__iter__()
    for i in range(num):
        it.__next__()


def compara_qname(read_name, mate_name):
    pattern = re.compile("^sq_[0-9]+_")
    read_match = re.search(pattern, read_name)
    mate_match = re.search(pattern, mate_name)
    read_sub = int(read_name[read_match.start()+3:read_match.end()-1])
    mate_sub = int(mate_name[mate_match.start()+3:mate_match.end()-1])
    if read_sub < mate_sub:
        return True
    else:
        return False


def compara_file():
    read_lette = 0
    bam_file = apri_bam_file(1)
    for read in bam_file:
        read_not_equal = 0
        mate_file = apri_bam_file(2)
        skip_read(mate_file, read_lette)
        for mate in mate_file:
            if getQueryName(read) == getQueryName(mate):
                insert_length.append(math.fabs(
                    read.reference_start - mate.reference_start))
                read_lette += 1
                if read_not_equal > 0:
                    read_lette += (read_not_equal)
                    print("read_not_equal ", read_not_equal,
                        "read_lette ", read_lette)
                print("equal")
                print(insert_length[len(insert_length) - 1])
                print(read.query_name, "\t", mate.query_name)
                break
            else:
                print("not equal")
                print(read.query_name, "\t", mate.query_name)
                if compara_qname(read.query_name, mate.query_name):
                    if read_not_equal > 0:
                        read_lette += (read_not_equal - 1)
                        print("if", "read_lette", read_lette,
                            "read_not_equal", read_not_equal)
                    break
                else:
                    read_not_equal += 1
                    print("else", read_not_equal)
        mate_file.close()
    bam_file.close()

if __name__ == "__main__":
    compara_file()
