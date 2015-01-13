#Resequencing project
#modulo principale

from reseq_stampa import *
from reseq_utility import *
import math

#massima lunghezza tollerata per la lunghezza degli inserti
max_insert_length = 20000


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


def esamina_bam_dict():
    """Esamina il file bam per ricercare le single, unique e multiple
    read presenti; una volta trovate, le stampa su file sam."""
    stampa_messaggi("inizio esamina")
    #apertura bam file
    bam_file = apri_bam_file("all")
    #dizionario dove mappare le read uguali
    read_dict = dict()
    #riferimenti alle liste relative alle single, unique e multiple read
    single_reads = unique_reads = multiple_reads = []
    for read in bam_file:
        name = get_query_name(read)
        if not name in read_dict:
            read_dict[name] = []
        read_dict[name].append(read)
        print(read_dict)
    for key, read_list in read_dict.items():
        if len(read_dict[key]) == 1:
            single_reads.append(read_list[0])
        elif len(read_dict[key]) == 2:
            unique_reads.append(read_list[0])
            unique_reads.append(read_list[1])
        else:
            for read in read_list:
                multiple_reads.append(read)
    stampa_messaggi("fine esamina")
    stampa_messaggi("inizio stampa file")
    stampa_read(single_reads, bam_file, "single_reads_dict.sam")
    stampa_read(unique_reads, bam_file, "unique_reads_dict.sam")
    stampa_read(multiple_reads, bam_file, "multiple_reads_dict.sam")
    bam_file.close()
    stampa_messaggi("fine stampa file")


def calcola_insert_length():
    pass


if __name__ == "__main__":
    #esamina_bam()
    esamina_bam_dict()
    #calcola_insert_length()
