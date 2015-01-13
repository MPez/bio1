#Resequencing project
#modulo con funzioni di utilità

import pysam
import re

bam_file_name1 = "pass_bam/pass_reads1_sorted_name.bam"
bam_file_name2 = "pass_bam/pass_reads2_sorted_name.bam"
bam_file_all = "pass_bam/pass_reads_all_sorted_name.bam"


def apri_bam_file(name):
    """Apre il file bam desiderato in base all'indice passato come argomento.

    Ritorna l'oggetto che corrisponde al file bam."""
    if name == 1:
        return pysam.AlignmentFile(bam_file_name1, "rb")
    elif name == 2:
        return pysam.AlignmentFile(bam_file_name2, "rb")
    elif name == "all":
        return pysam.AlignmentFile(bam_file_all, "rb")


def get_query_name(query):
    """Calcola la query name della read passata come argomento.

    Ritorna la stringa al netto dei caratteri finali che identificano
    il file di provenienza."""
    qname_pattern = re.compile("/[1|2]$")
    return re.sub(qname_pattern, "", query.query_name)
