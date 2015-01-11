#Resequencing project

import pysam
import re
import math
import statistics
from multiprocessing import Pool

bam_file_name1 = "pass_bam/pass_reads1_sorted_name.bam"
bam_file_name2 = "pass_bam/pass_reads2_sorted_name.bam"
max_insert_length = 20000


def apri_bam_file(num):
    """Apre il file bam desiderato in base all'indice passato come argomento.

    Ritorna l'oggetto che corrisponde al file bam."""
    if num == 1:
        return pysam.AlignmentFile(bam_file_name1, "rb")
    else:
        return pysam.AlignmentFile(bam_file_name2, "rb")


def get_query_name(query):
    """Calcola la query name della read passata come argomento.

    Ritorna la stringa al netto dei caratteri finali che identificano
    il file di provenienza."""
    qname_pattern = re.compile("/[1|2]$")
    return re.sub(qname_pattern, "", query.query_name)


def compara_qname(read_name, mate_name):
    """Compara le query name delle 2 read passate come argomento;
    ogni query name comprende le posizioni di entrambe le read:
    se sono uguali corrispondono a 2 mate pair, altrimenti si deve
    controllare quale read è singola.

    Ritorna True nel caso read_name preceda mate_name;
    Ritorna False nel caso read_name segua mate_name."""
    #pattern per ottenere la posizione della prima read nelle 2 query
    pattern = re.compile("^sq_[0-9]+_")
    read_match = re.search(pattern, read_name)
    mate_match = re.search(pattern, mate_name)
    read_sub = int(read_name[read_match.start()+3:read_match.end()-1])
    mate_sub = int(mate_name[mate_match.start()+3:mate_match.end()-1])
    #se la poszione della prima read è minore nella prima query
    if read_sub < mate_sub:
        return True
    #se le posizioni della prima read coincide nelle 2 query
    elif read_sub == mate_sub:
        #pattern per ottenere la posizione della seconda read nelle 2 query
        pattern2 = re.compile("_[0-9]+_[0|1]_")
        read2_match = re.search(pattern2, read_name)
        mate2_match = re.search(pattern2, mate_name)
        read2_sub = int(read_name[read2_match.start()+1:read2_match.end()-3])
        mate2_sub = int(mate_name[mate2_match.start()+1:mate2_match.end()-3])
        #se la posizione della seconda read è minore nella prima query
        if read2_sub < mate2_sub:
            return True
    #se la posizione della prima o della seconda read
    #è minore nella seconda query
    else:
        return False


def stampa_length(insert_list, discarded_list):
    """Crea e scrive file di testo relativi ai mate pair trovati e
    alla lunghezza degli inserti corrispondenti.
    Tali file serviranno per essere usati con gnuplot."""
    insert_file = open("insert.txt", "w")
    discarded_insert_file = open("discarded_insert.txt", "w")
    val = 0
    for (i, j) in insert_list:
        val += 1
        insert_file.write(str(val) + "\t" + str(i) + "\t" + str(j) + "\n")
    val = 0
    for (i, j) in discarded_list:
        val += 1
        discarded_insert_file.write(str(val) + "\t" +
                                    str(i) + "\t" + str(j) + "\n")
    insert_file.close()
    discarded_insert_file.close()


def stampa_unique(unique_read, bam_file):
    """Crea e stampa un sam file dove vengono scritte le unique read trovate.
    """
    unique_file = pysam.AlignmentFile("unique_reads.sam", "w",
                                      referencenames=bam_file.references,
                                      referencelengths=bam_file.lengths)
    for read in unique_read:
        unique_file.write(read)
    unique_file.close()


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


def calcola_isert_length():
    """Compara le read dei 2 file sam in esame trovando i mate pair:
    i mate pair posseggono la stessa query name a meno dei 2 caratteri finali
    che identificano il file bam di provenienza."""
    stampa_messaggi("inizio lenght")
    #apertura file bam
    bam_file = apri_bam_file(1)
    mate_file = apri_bam_file(2)
    #creazione iteratori
    bam_it = iter(bam_file)
    mate_it = iter(mate_file)
    #riferimenti alle prime read dei file
    read = next(bam_it)
    mate = next(mate_it)
    #riferimento alla read precedente
    #previous_read = None
    #liste di tuple contenenti nome read e lunghezza inserti
    insert_length = []
    discarded_insert_length = []
    #liste contenenti le read uniche e multiple trovate
    unique_read = []
    #multiple_read = []
    #variabile di controllo ciclo
    terminato = False
    while not terminato:
        try:
            #se le query name coincidono le read sono mate pair
            if get_query_name(read) == get_query_name(mate):
                length = math.fabs(read.reference_start - mate.reference_start)
                #vengono scartati gli inserti con lunghezze fuori range
                if length < max_insert_length:
                    insert_length.append((get_query_name(read), length))
                else:
                    discarded_insert_length.append(
                        (get_query_name(read), length))
                #viene aggiornata la read precedente che ha trovato il mate
                #previous_read = read
                #avanzano entrambi gli iteratori dei file bam
                read = next(bam_it)
                mate = next(mate_it)
            #se la prima query precede la seconda,
            #è singola e deve essere saltata
            elif compara_qname(read.query_name, mate.query_name):
                unique_read.append(read)
                read = next(bam_it)
            #se la seconda query precede la prima,
            #è singola e deve essere saltata
            else:
                unique_read.append(mate)
                mate = next(mate_it)
        except StopIteration:
            terminato = True
    stampa_messaggi("fine lenght")
    stampa_length(insert_length, discarded_insert_length)
    stampa_unique(unique_read, bam_file)
    stampa_messaggi("stat", insert_length, discarded_insert_length)
    bam_file.close()
    mate_file.close()


def stampa_multiple(multiple_read):
    multiple_file = open("multiple_read.txt", "w")
    for (seq, name) in multiple_read:
        multiple_file.write(seq + "\t" + name + "\n")
    multiple_file.close()


def cerca_query(read, num):
    multiple_read = []
    seq = read.query_sequence
    bam_file1 = apri_bam_file(num)
    for read1 in bam_file1:
        if seq == read1.query_sequence and read.query_name != read1.query_name:
            multiple_read.append((seq, read1.query_name))
    bam_file1.close()
    return multiple_read


def conta_multiple_read(num):
    stampa_messaggi("inizio multiple")
    multiple_read = []
    bam_file = apri_bam_file(num)
    pass
    bam_file.close()
    if len(multiple_read) > 0:
        stampa_multiple(multiple_read)
    else:
        print("non sono state trovate read multiple")
    stampa_messaggi("fine multiple")


if __name__ == "__main__":
    calcola_isert_length()
    #conta_multiple_read(1)
    #conta_multiple_read(2)
