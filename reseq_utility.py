# Resequencing project
# Pezzutti Marco - 1084411
# modulo con funzioni di utilità

import pysam
import re
import os
from reseq_stampa import dir_risultati


# cartella di destinazione dei file pass bam forniti
bam_dir = None
# riferimenti ai nomi dei file bam
bam_file_name1 = "pass_reads1_sorted_name.bam"
bam_file_name2 = "pass_reads2_sorted_name.bam"
bam_file_all = "pass_reads_all_sorted_name.bam"
unique_file = "unique_reads.bam"
single_file = "single_reads.bam"
multiple_file = "multiple_reads.bam"
unique_sorted = "unique_reads_sorted.bam"
single_sorted = "single_reads_sorted.bam"
multiple_sorted = "multiple_reads_sorted.bam"


def set_bam_dir(dir):
    """Imposta la variabile globale relativa alla cartella
    di destinazione dei file bam."""
    global bam_dir
    bam_dir = dir


def chiama_reseq():
    """Richiama lo script bash reseq ed esegue il sort dei file bam creati.
    """
    os.system("bash reseq.sh -t")


def apri_bam_file(name):
    """Apre il file bam desiderato in base all'indice passato come argomento.

    Ritorna l'oggetto che corrisponde al file bam."""
    if name == 1:
        return pysam.AlignmentFile(bam_dir + bam_file_name1, "rb")
    elif name == 2:
        return pysam.AlignmentFile(bam_dir + bam_file_name2, "rb")
    elif name == "all":
        return pysam.AlignmentFile(bam_dir + bam_file_all, "rb")
    elif name == "single":
        return pysam.AlignmentFile(dir_risultati + single_file, "rb")
    elif name == "unique":
        return pysam.AlignmentFile(dir_risultati + unique_file, "rb")
    elif name == "multiple":
        return pysam.AlignmentFile(dir_risultati + multiple_file, "rb")
    elif name == "single_sorted":
        return pysam.AlignmentFile(dir_risultati + single_sorted, "rb")
    elif name == "unique_sorted":
        return pysam.AlignmentFile(dir_risultati + unique_sorted, "rb")
    elif name == "multiple_sorted":
        return pysam.AlignmentFile(dir_risultati + multiple_sorted, "rb")


def get_query_name(query):
    """Calcola la query name della read passata come argomento.

    Ritorna la stringa al netto dei caratteri finali che identificano
    il file di provenienza."""
    qname_pattern = re.compile("/[1|2]$")
    return re.sub(qname_pattern, "", query.query_name)


def find_number(read):
    """Cerca se la read passata come argomento appartiene al primo
    file pass o al secondo fornito.

    Ritorna true se appartiene al secondo e false al primo."""
    one = re.compile("/1$")
    if re.search(one, read.query_name):
        return True
    else:
        return False
