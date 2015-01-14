#Resequencing project
#modulo con funzioni di stampa

import pysam
import statistics


def stampa_messaggi(mess, insert=0, discarded=0):
    """Stampa messaggi di testo per l'utente con informazioni
    sull'avanzamento e sulle statistiche ottenute."""
    if mess == "inizio lenght":
        print("calcolo lunghezza inserti iniziato")
    elif mess == "fine length":
        print("calcolo lunghezza inserti terminato")
        print("preparazione file input per gnuplot")
    elif mess == "stat":
        print("sono stati rilevati", len(insert), "mate pair")
        print("sono stati scartati", len(discarded),
              "mate pair in quanto fuori range")
        print("valori statistici sulla lunghezza degli inserti genomici")
        print("media", statistics.mean(j for (i, j) in insert),
              "mediana", statistics.median(j for (i, j) in insert),
              "varianza", statistics.variance(j for (i, j) in insert),
              "deviazione standard", statistics.stdev(j for (i, j) in insert))
    elif mess == "inizio esamina":
        print("analisi bam file per ricerca single, "
              "unique e multiple reads iniziata")
    elif mess == "fine esamina":
        print("analisi bam file terminata")
    elif mess == "inizio stampa file":
        print("creazione bam file con single, "
              "unique e multiple reads iniziata")
    elif mess == "fine stampa file":
        print("creazione bam file terminata")


def stampa_read(read_list, bam_file, nome_file):
    """Crea e stampa un sam file dove vengono scritte le read trovate.
    """
    new_file = pysam.AlignmentFile(nome_file + ".bam", "wb",
                                   header=bam_file.header)
    for read in read_list:
        new_file.write(read)
    new_file.close()


def stampa_length(read_list, nome_file):
    """Crea e scrive un file di testo relativo ai mate pair trovati e
    alla lunghezza degli inserti corrispondenti.
    Tale file serve per essere usati con gnuplot."""
    read_file = open(nome_file, "w")
    val = 0
    for (nome, length) in read_list:
        val += 1
        read_file.write(str(val) + "\t" + str(nome) + "\t"
                        + str(length) + "\n")
    read_file.close()


def stampa_coverage(coverage, file_name, region_name):
    new_file = open(file_name + ".wig", "w")
    new_file.write("variableStep chrom=" + region_name + "\n")
    for (pos, n) in coverage:
        new_file.write(str(pos) + "\t" + str(n) + "\n")
    new_file.close()
