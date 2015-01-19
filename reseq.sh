#!/bin/bash
# Resequencing project
# Pezzutti Marco - 1084411

ris_dir="risultati"
bam_dir="pass_bam"

# avvia lo script python che esegue il resequencing
start_reseq() {
    script="reseq_new.py"
    python $script
}


# converte tutti i file bam presenti nella cartella e 
# li converte in file sam
bam_to_sam() {
    echo "conversione da file bam a file sam"
    suff="bam"
    SUFF="sam"
    cd $1
    for s in $(ls *.$suff); do
        samtools view -o ${s%.$suff}.$SUFF -h $s
    done
    cd ..
}


# effettua il merge dei file pass, li ordina per nome e per posizione e crea l'indice
merge_pass_sort() {
    suff=".bam"
    pass1="pass_reads1"
    pass2="pass_reads2"
    all_name="pass_reads_all"
    all_sorted="pass_reads_all_sorted"
    name="_name"
    temp="temp"

    echo "merge dei file pass $pass1$suff e $pass2$suff"
    cd $bam_dir
    samtools merge $all_name$suff $pass1$suff $pass2$suff

    echo "sort per nome read del file $all_name$suff"
    samtools sort -o $all_sorted$name$suff -n -T $temp$suff $all_name$suff

    echo "sort per posizione read e creazione indice del file $all_name$suff"
    samtools sort -o $all_sorted$suff -T $temp$suff $all_name$suff
    samtools index -b $all_sorted$suff
    cd ..
}


# crea dai file wig i relativi file tdf
wig_to_tdf() {
    echo "creazione file tdf"
    suff="wig"
    SUFF="tdf"
    genome="./reference_A_laidlawii/reference_A_laidlawii.fasta.fai"

    cd $ris_dir
    for s in $(ls *.$suff); do
        igvtools toTDF $s $s.$SUFF $genome
    done
    cd ..
}


# genera i grafici con gnuplot
plot() {
    echo "generazione grafici lunghezza inserti"
    graph="plot.gnuplot"
    gnuplot $graph
}

# main
while getopts ":a:b:gmsw" flag; do
    case $flag in
        a )
            merge_pass_sort
            start_reseq
            bam_to_sam $OPTARG
            wig_to_tdf
            gnuplot
            ;;
        b ) bam_to_sam $OPTARG
            ;;
        g ) plot
            ;;
        m ) merge_pass_sort
            ;;
        s ) start_reseq
            ;;
        w ) wig_to_tdf
            ;;
        * ) echo "Opzione non presente, nessuna operazione da svolgere."
            ;;
    esac
done