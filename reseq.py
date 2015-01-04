#Resequencing project

import pysam
import re
import math

bam_file_name = "pass_bam/pass_reads1_sorted_name.bam"
mate_file_name = "pass_bam/pass_reads2_sorted_name.bam"
pattern = re.compile("/[1|2]$")
insert_length = []


def apri_bam_file(name):
    if name == "bam":
        return pysam.AlignmentFile(bam_file_name, "rb")
    else:
        return pysam.AlignmentFile(mate_file_name, "rb")


def get_iterator(iterable):
    return iterable.__iter__()


def prossima_read(iterator):
    return iterator.__next__()


def get_query_name(query):
    return re.sub(pattern, "", query.query_name)


def calcola_insert_length():
    bam_file = apri_bam_file("bam")
    mate_file = apri_bam_file("mate")
    read_it = get_iterator(bam_file)
    mate_it = get_iterator(mate_file)
    read = prossima_read(read_it)
    mate = prossima_read(mate_it)
    alt = 0
    terminato = False
    while not terminato:
        try:
            read_length = read.reference_start + 1
            mate_length = mate.reference_start + 1
            if get_query_name(read) == get_query_name(mate):
                insert_length.append(math.fabs(read_length - mate_length))
                print("equal")
                print(read.query_name, "\t", mate.query_name)
                print(insert_length[len(insert_length) - 1])
                read = prossima_read(read_it)
                mate = prossima_read(mate_it)
            else:
                if alt % 2:
                    read = prossima_read(read_it)
                else:
                    mate = prossima_read(mate_it)
                alt += 1
        except StopIteration:
            terminato = True
    bam_file.close()
    mate_file.close()


if __name__ == "__main__":
    calcola_insert_length()
