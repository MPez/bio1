#Resequencing project
#Pezzutti Marco - 1084411
#modulo principale

from reseq_stampa import *
from reseq_utility import *

#massima lunghezza tollerata per la lunghezza degli inserti
max_insert_length = 20000


def esamina_multiple(multiple_list, mate_list, single_list):
    for i in range(0, len(multiple_list)):
        trovato = False
        read = multiple_list[i]
        for j in range(i+1, len(multiple_list)):
            if not trovato:
                mate = multiple_list[j]
                if read.is_reverse and not mate.is_reverse:
                    if find_number(read) and not find_number(mate):
                        mate_list.append(read)
                        mate_list.append(mate)
                        trovato = True
                    elif not find_number(read) and find_number(mate):
                        mate_list.append(read)
                        mate_list.append(mate)
                        trovato = True
                elif not read.is_reverse and mate.is_reverse:
                    if find_number(read) and not find_number(mate):
                        mate_list.append(read)
                        mate_list.append(mate)
                        trovato = True
                    elif not find_number(read) and find_number(mate):
                        mate_list.append(read)
                        mate_list.append(mate)
                        trovato = True
        if j == len(multiple_list) - 1 and read not in mate_list:
            single_list.append(read)


def esamina_bam():
    """Esamina il file bam per ricercare le single, unique e multiple
    read presenti; una volta trovate, le stampa su file bam."""
    stampa_messaggi("inizio esamina")
    #apertura bam file
    bam_file = apri_bam_file("all")
    #creazione iteratore
    bam_it = iter(bam_file)
    #riferimenti alla read precedente e a quella corrente
    prev_read = None
    read = next(bam_it)
    #riferimenti alle liste relative alle single, unique e multiple read
    single_reads = []
    unique_reads = []
    multiple_reads = []
    multiple_mate_reads = []
    multiple_single_reads = []
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
                    equal_reads.append(prev_read)
                    for r in equal_reads:
                        multiple_reads.append(r)
                    esamina_multiple(equal_reads,
                                     multiple_mate_reads,
                                     multiple_single_reads)
                    equal_reads.clear()
                    count = 0
        except StopIteration:
            terminato = True
            #all'uscita dal ciclo rimangono da esaminare le ultime 2 read
            #presenti nel file bam di origine, vengono esaminate come le altre
            if len(equal_reads) > 0:
                if count == 0:
                    single_reads.append(prev_read)
                elif count == 1:
                    unique_reads.append(equal_reads.pop())
                    unique_reads.append(prev_read)
                    count = 0
                else:
                    equal_reads.append(prev_read)
                    for r in equal_reads:
                        multiple_reads.append(r)
                    esamina_multiple(equal_reads,
                                     multiple_mate_reads,
                                     multiple_single_reads)
                    equal_reads.clear()
    #serie di comandi di stampa informazioni per l'utente e file bam
    stampa_messaggi("fine esamina")
    stampa_messaggi("inizio stampa file")
    stampa_read(single_reads, bam_file, "single_reads")
    stampa_read(unique_reads, bam_file, "unique_reads")
    stampa_read(multiple_reads, bam_file, "multiple_reads")
    stampa_read(multiple_mate_reads, bam_file, "multiple_mate_reads")
    stampa_read(multiple_single_reads, bam_file, "multiple_single_reads")
    bam_file.close()
    stampa_messaggi("fine stampa file")


def calcola_insert_length():
    """Calcola la lunghezza degli inserti composti dagli unique mate pair;
    li filtra se sono fuori range massimo e li stampa su due file separati.
    """
    stampa_messaggi("inizio length")
    #apertura bam file
    unique_file = apri_bam_file("unique")
    #creazione iteratore
    unique_it = iter(unique_file)
    #riferimenti alla read e al suo mate
    read = next(unique_it)
    mate = next(unique_it)
    #riferimenti a caratteristiche del genoma in esame
    region_number = unique_file.nreferences
    region_name = unique_file.references[region_number - 1]
    region_length = unique_file.lengths[region_number - 1]
    #riferimento alla lista relativa alla physical coverage
    coverage = [[pos, 0] for pos in range(1, region_length + 1)]
    #riferimenti alle liste relative a lunghezza degli inserti
    insert_length = []
    discarded_insert_length = []
    insert_len_wig = [[pos, 0] for pos in range(1, region_length + 1)]
    disc_insert_len_wig = [[pos, 0] for pos in range(1, region_length + 1)]
    #variabile di controllo del ciclo
    terminato = False
    while not terminato:
        try:
            #posizione di inizio delle read
            read_pos = read.reference_start
            mate_pos = mate.reference_start
            #controllo sulla posizione delle read
            if read_pos < mate_pos:
                length = mate_pos + mate.reference_length - read_pos
                #se la lunghezza dell'inserto è nel range,
                #inserisco i dati nelle liste
                if length < max_insert_length:
                    insert_len_wig[read_pos][1] += length
                    insert_len_wig[mate_pos][1] += length
                    insert_length.append((get_query_name(read), length))
                    for i in range(read_pos, length):
                        coverage[i][1] += 1
                #se la lunghezza è fuori range escludo i mate pair
                else:
                    disc_insert_len_wig[read_pos][1] += length
                    disc_insert_len_wig[mate_pos][1] += length
                    discarded_insert_length.append((get_query_name(read),
                                                    length))
            #procedimento analogo al precedente con ordine delle read invertite
            else:
                length = read_pos + read.reference_length - mate_pos
                if length < max_insert_length:
                    insert_len_wig[mate_pos][1] += length
                    insert_len_wig[read_pos][1] += length
                    insert_length.append((get_query_name(read), length))
                    for i in range(mate_pos, read_pos + read.reference_length):
                        coverage[i][1] += 1
                else:
                    disc_insert_len_wig[mate_pos][1] += length
                    disc_insert_len_wig[read_pos][1] += length
                    discarded_insert_length.append((get_query_name(read),
                                                    length))
            #avanzo gli iteratori per ottenere i nuovi mate pair da controllare
            read = next(unique_it)
            mate = next(unique_it)
        except StopIteration:
            terminato = True
    #serie di messaggi informativi e di stampa liste in file gnuplot e wiggle
    stampa_messaggi("fine length")
    stampa_messaggi("inizio stampa gnuplot")
    stampa_length(insert_length, "insert_length.txt")
    stampa_messaggi("fine stampa gnuplot")
    stampa_length(discarded_insert_length, "discarded_insert_length.txt")
    stampa_messaggi("inizio stampa wiggle")
    stampa_coverage(coverage, "unique_physical", region_name)
    stampa_coverage(insert_len_wig, "insert_length", region_name)
    stampa_coverage(disc_insert_len_wig, "discarded_insert_length",
                    region_name)
    stampa_messaggi("fine stampa wiggle")
    stampa_messaggi("stat", insert_length, discarded_insert_length)
    unique_file.close()


def calcola_coverage(read_type):
    """Calcola la sequence coverage del tipo di read richiesto e la stampa
    su un file wiggle."""
    stampa_messaggi("inizio coverage", read_type=read_type)
    #apertura bam file
    bam_file = apri_bam_file(read_type)
    #riferimenti a caratteristiche del genoma in esame
    region_number = bam_file.nreferences
    region_name = bam_file.references[region_number - 1]
    region_length = bam_file.lengths[region_number - 1]
    #refiremento alla lista relativa alla sequence coverage
    coverage = []
    #per ogni base del genoma viene calcolata la coverage
    for col in bam_file.pileup(region_name, 0, region_length - 1):
        coverage.append((col.reference_pos, col.nsegments))
    #serie di messaggi informativi e di stampa lista in file wiggle
    stampa_messaggi("fine coverage")
    stampa_messaggi("inizio stampa wiggle")
    stampa_coverage(coverage, read_type, region_name)
    stampa_messaggi("fine stampa wiggle")
    bam_file.close()


if __name__ == "__main__":
    #esamina_bam()
    calcola_insert_length()
    #calcola_coverage("unique_sorted")
    #calcola_coverage("single_sorted")
    #calcola_coverage("multiple_sorted")
