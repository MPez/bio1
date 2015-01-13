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
    elif mess == "inizio multiple":
        print("calcolo read multiple iniziato")
    elif mess == "fine multiple":
        print("calcolo read multiple terminato")
    elif mess == "inizio esamina":
        print("analisi bam file per ricerca single, unique e multiple\
              reads iniziata")
    elif mess == "fine esamina":
        print("analisi bam file terminata")
    elif mess == "inizio stampa file":
        print("creazione sam file con signle, unique e multiple reads\
              iniziata")
    elif mess == "fine stampa file":
        print("creazione sam file terminata")


def stampa_read(read_list, bam_file, nome_file):
    """Crea e stampa un sam file dove vengono scritte le read trovate.
    """
    new_file = pysam.AlignmentFile(nome_file, "w",
                                   referencenames=bam_file.references,
                                   referencelengths=bam_file.lengths)
    for read in read_list:
        new_file.write(read)
    new_file.close()
