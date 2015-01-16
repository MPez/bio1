#!/bin/bash
#Resequencing project
#Pezzutti Marco - 1084411


#avvia lo script python che esegue il resequencing
script="reseq_new.py"
python $script


#converte tutti i file bam presenti nella cartella e 
#li converte in file sam
if [ $1 =="sam" ]; then
    echo "conversione da file bam a file sam"
    suff="bam"
    SUFF="sam"
    for s in $(ls *.$suff); do
        samtools view -o ${s%.$suff}.$SUFF -h $s
    done
fi


#crea dai file wig i relativi file tdf
echo "creazione file tdf"
suff="wig"
SUFF="tdf"
genome="./reference_A_laidlawii/reference_A_laidlawii.fasta.fai"
for s in $(ls *.$suff); do
    igvtools toTDF $s $s.$SUFF $genome
done


#genera i grafici con gnuplot
echo "generazione grafici lunghezza inserti"
graph="plot.gnuplot"
gnuplot $graph
