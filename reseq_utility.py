#Resequencing project
#modulo con funzioni di utilità

import pysam
import re

#riferimenti ai nomi dei file bam
bam_file_name1 = "pass_bam/pass_reads1_sorted_name.bam"
bam_file_name2 = "pass_bam/pass_reads2_sorted_name.bam"
bam_file_all = "pass_bam/pass_reads_all_sorted_name.bam"
unique_file = "unique_reads.bam"
single_file = "single_reads.bam"
multiple_file = "multiple_reads.bam"
unique_sorted = "unique_reads_sorted.bam"
single_sorted = "single_reads_sorted.bam"
multiple_sorted = "multiple_reads_sorted.bam"


def apri_bam_file(name):
    """Apre il file bam desiderato in base all'indice passato come argomento.

    Ritorna l'oggetto che corrisponde al file bam."""
    if name == 1:
        return pysam.AlignmentFile(bam_file_name1, "rb")
    elif name == 2:
        return pysam.AlignmentFile(bam_file_name2, "rb")
    elif name == "all":
        return pysam.AlignmentFile(bam_file_all, "rb")
    elif name == "single":
        return pysam.AlignmentFile(single_file, "rb")
    elif name == "unique":
        return pysam.AlignmentFile(unique_file, "rb")
    elif name == "multiple":
        return pysam.AlignmentFile(multiple_file, "rb")
    elif name == "single_sorted":
        return pysam.AlignmentFile(single_sorted, "rb")
    elif name == "unique_sorted":
        return pysam.AlignmentFile(unique_sorted, "rb")
    elif name == "multiple_sorted":
        return pysam.AlignmentFile(multiple_sorted, "rb")


def get_query_name(query):
    """Calcola la query name della read passata come argomento.

    Ritorna la stringa al netto dei caratteri finali che identificano
    il file di provenienza."""
    qname_pattern = re.compile("/[1|2]$")
    return re.sub(qname_pattern, "", query.query_name)
