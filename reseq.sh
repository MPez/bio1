#!/bin/bash
#Resequencing project
#Pezzutti Marco - 1084411


#avvia lo script python che esegue il resequencing
start_reseq() {
    script="reseq_new.py"
    python $script
}


#converte tutti i file bam presenti nella cartella e 
#li converte in file sam
bam_to_sam() {
    echo "conversione da file bam a file sam"
    suff="bam"
    SUFF="sam"
    for s in $(ls *.$suff); do
        samtools view -o ${s%.$suff}.$SUFF -h $s
    done
}


#crea dai file wig i relativi file tdf
wig_to_tdf() {
    echo "creazione file tdf"
    suff="wig"
    SUFF="tdf"
    genome="./reference_A_laidlawii/reference_A_laidlawii.fasta.fai"
    for s in $(ls *.$suff); do
        igvtools toTDF $s $s.$SUFF $genome
    done
}


#genera i grafici con gnuplot
gnuplot() {
    echo "generazione grafici lunghezza inserti"
    graph="plot.gnuplot"
    gnuplot $graph
}

#main
while getopts ":abgsw" flag; do
    case $flag in
        a )
            start_reseq
            bam_to_sam
            wig_to_tdf
            gnuplot
            ;;
        b ) bam_to_sam
            ;;
        g ) gnuplot
            ;;
        s ) start_reseq
            ;;
        w ) wig_to_tdf
            ;;
        * ) echo "Opzione non presente, nessuna operazione da svolgere."
            ;;
    esac
done