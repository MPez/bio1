#Resequencing project

import pysam
import re
import math

bamFileName = "pass_bam/pass_reads_all_sorted_name.bam"
pattern = re.compile("/[1|2]$")
insertLength = []


def apriBamFile():
    return pysam.AlignmentFile(bamFileName, "rb")


def getIterator(iterable):
    return iterable.__iter__()


def prossimaRead(iterator):
    return iterator.__next__()


def prossimaCoppia(iterator):
    prossimaRead(iterator)
    return prossimaRead(iterator)


def getQueryName(query):
    return re.sub(pattern, "", query.query_name)


def calcolaInsertLength():
    bamFile = apriBamFile()
    mateFile = apriBamFile()
    readIt = getIterator(bamFile)
    mateIt = getIterator(mateFile)
    alt = 0
    terminato = False
    while not terminato:
        try:
            read = prossimaRead(readIt)
            readLength = read.reference_start + 1
            mate = prossimaRead(mateIt)
            mateLength = mate.reference_start + 1
            if getQueryName(read) == getQueryName(mate):
                insertLength.append(math.fabs(readLength - mateLength))
                print("equal")
                print(read.query_name, "\t", mate.query_name)
                print(insertLength[len(insertLength) - 1])
                read = prossimaCoppia(readIt)
                mate = prossimaCoppia(mateIt)
            else:
                if alt % 2:
                    read = prossimaRead(readIt)
                else:
                    mate = prossimaRead(mateIt)
        except StopIteration:
            terminato = True
    bamFile.close()
    mateFile.close()


if __name__ == "__main__":
    calcolaInsertLength()
