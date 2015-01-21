# Resequencing project
# Pezzutti Marco - 1084411
# modulo con funzioni di stampa

import pysam
import statistics
import re
import os

# cartella dove vengono salvati i file prodotti
global dir_risultati
dir_risultati = "risultati/"

# crea la cartella di destinazione dei risultati se non esiste
if not os.path.exists(dir_risultati):
    os.makedirs(dir_risultati)


def stampa_messaggi(mess, insert=0, discarded=0, read_type=0):
    """Stampa messaggi di testo per l'utente con informazioni
    sull'avanzamento e sulle statistiche ottenute."""
    if mess == "inizio length":
        print("calcolo lunghezza inserti e physical coverage iniziato")
    elif mess == "fine length":
        print("calcolo lunghezza inserti e physical coverage terminato")
    elif mess == "inizio stampa gnuplot":
        print("preparazione file input per gnuplot")
    elif mess == "fine stampa gnuplot":
        print("preparazione file input per gnuplot terminata")
    elif mess == "stat":
        print("sono stati rilevati", len(insert), "unique pair")
        print("sono stati scartati", len(discarded),
              "mate pair in quanto fuori range")
        print("valori statistici sulla lunghezza degli inserti genomici")
        print("media", statistics.mean(j for (i, j) in insert),
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
    elif mess == "inizio coverage":
        name = re.sub(re.compile("_[a-z]+$"), "", read_type)
        print("calcolo sequence coverage per %s read" % (name))
    elif mess == "fine coverage":
        print("calcolo sequence coverage termianto")
    elif mess == "inizio stampa wiggle":
        print("creazione file wiggle iniziata")
    elif mess == "fine stampa wiggle":
        print("creazione file wiggle terminata")


def stampa_read(read_list, bam_file, nome_file):
    """Crea e stampa un sam file dove vengono scritte le read trovate.
    """
    new_file = pysam.AlignmentFile(dir_risultati + nome_file + ".bam", "wb",
                                   header=bam_file.header)
    for read in read_list:
        new_file.write(read)
    new_file.close()


def stampa_length(read_list, nome_file):
    """Crea e scrive un file di testo relativo ai mate pair trovati e
    alla lunghezza degli inserti corrispondenti.
    Tale file serve per essere usati con gnuplot."""
    read_file = open(dir_risultati + nome_file, "w")
    val = 0
    for (nome, length) in read_list:
        val += 1
        read_file.write(str(val) + "\t" + str(nome) + "\t"
                        + str(length) + "\n")
    read_file.close()


def stampa_coverage(coverage, file_name, region_name):
    """Crea e scrive un file wiggle relativo alla coverage calcolata.
    """
    file_name = re.sub(re.compile("_[a-z]+$"), "", file_name)
    new_file = open(dir_risultati + file_name + "_coverage.wig", "w")
    new_file.write("variableStep chrom=" + region_name + "\n")
    for (pos, n) in coverage:
        new_file.write(str(pos) + "\t" + str(n) + "\n")
    new_file.close()
