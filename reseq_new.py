#Resequencing project

import pysam
import re
import math
import statistics

bam_file_name1 = "pass_bam/pass_reads1_sorted_name.bam"
bam_file_name2 = "pass_bam/pass_reads2_sorted_name.bam"
bam_file_all = "pass_bam/pass_reads_all_sorted_name.bam"
max_insert_length = 20000


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


def stampa_single(single_read, bam_file):
    """Crea e stampa un sam file dove vengono scritte le unique read trovate.
    """
    single_file = pysam.AlignmentFile("single_reads.sam", "w",
                                      referencenames=bam_file.references,
                                      referencelengths=bam_file.lengths)
    for read in single_read:
        single_file.write(read)
    single_file.close()


def stampa_read(read_list, bam_file, nome_file):
    """Crea e stampa un sam file dove vengono scritte le read trovate.
    """
    new_file = pysam.AlignmentFile(nome_file, "w",
                                   referencenames=bam_file.references,
                                   referencelengths=bam_file.lengths)
    for read in read_list:
        new_file.write(read)
    new_file.close()


def esamina_bam():
    """Esamina il file bam per ricercare le single, unique e multiple
    read presenti; una volta trovate, le stampa su file sam."""
    stampa_messaggi("inizio esamina")
    #apertura bam file
    bam_file = apri_bam_file("all")
    #creazione iteratore
    bam_it = iter(bam_file)
    #riferimenti alla read precedente e a quella corrente
    prev_read = None
    read = next(bam_it)
    #riferimenti alle liste relative alle single, unique e multiple read
    single_reads = unique_reads = multiple_reads = []
    #riferimento alla lista dove vengono salvate le read uguali trovate
    equal_reads = []
    #variabile di controllo del ciclo
    terminato = False
    #contatore delle read uguali trovate
    count = 0
    while not terminato:
        try:
            #nuove read da controllare
            prev_read = read
            read = next(bam_it)
            #se le read sono uguali aggiorno il contatore e aggiungo alla lista
            if get_query_name(prev_read) == get_query_name(read):
                count += 1
                equal_reads.append(prev_read)
            #se le read non sono uguali controllo il contatore e salvo
            #le read trovate nella lista corretta
            else:
                if count == 0:
                    single_reads.append(prev_read)
                elif count == 1:
                    unique_reads.append(equal_reads.pop())
                    unique_reads.append(prev_read)
                    count = 0
                else:
                    for r in equal_reads:
                        multiple_reads.append(r)
                    multiple_reads.append(prev_read)
                    count = 0
                    equal_reads.clear()
        except StopIteration:
            terminato = True
    stampa_messaggi("fine esamina")
    stampa_messaggi("inizio stampa file")
    stampa_read(single_reads, bam_file, "single_reads.sam")
    stampa_read(unique_reads, bam_file, "unique_reads.sam")
    stampa_read(multiple_reads, bam_file, "multiple_reads.sam")
    bam_file.close()
    stampa_messaggi("fine stampa file")

if __name__ == "__main__":
    #calcola_insert_length()
    #conta_multiple_read("all")
    esamina_bam()
