#Resequencing project

import pysam
import re
import math
import statistics

bam_file_name1 = "pass_bam/pass_reads1_sorted_name.bam"
bam_file_name2 = "pass_bam/pass_reads2_sorted_name.bam"


def apri_bam_file(num):
    if num == 1:
        return pysam.AlignmentFile(bam_file_name1, "rb")
    else:
        return pysam.AlignmentFile(bam_file_name2, "rb")


def get_query_name(query):
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


def get_iterator(file):
    return file.__iter__()


def next_read(iterator):
    return iterator.__next__()


def compara():
    print("comparazione bam file iniziata")
    bam_file = apri_bam_file(1)
    mate_file = apri_bam_file(2)
    bam_it = get_iterator(bam_file)
    mate_it = get_iterator(mate_file)
    read = next_read(bam_it)
    mate = next_read(mate_it)
    insert_length = []
    terminato = False
    while not terminato:
        try:
            if get_query_name(read) == get_query_name(mate):
                insert_length.append(math.fabs(
                    read.reference_start - mate.reference_start))
                read = next_read(bam_it)
                mate = next_read(mate_it)
            else:
                if compara_qname(read.query_name, mate.query_name):
                    read = next_read(bam_it)
                else:
                    mate = next_read(mate_it)
        except StopIteration:
            terminato = True
    print("comparazione bam file completata")
    print("sono stati rilevati", len(insert_length) - 1, "mate pair")
    print("valori statistici sulla lunghezza degli inserti genomici")
    print("media", statistics.mean(insert_length),
          "mediana", statistics.median(insert_length),
          "varianza", statistics.variance(insert_length),
          "deviazione standard", statistics.stdev(insert_length))

if __name__ == "__main__":
    compara()
