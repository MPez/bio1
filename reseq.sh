#!/bin/bash
# Resequencing project
# Pezzutti Marco - 1084411

ris_dir="risultati/"
bam_dir="pass_bam/"

# avvia lo script python che esegue il resequencing
start_reseq() {
    script="reseq_new.py"
    python $script $1
}


# converte tutti i file bam presenti nella cartella e 
# li converte in file sam
bam_to_sam() {
    echo "Conversione da file bam a file sam"
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

    echo "Merge dei file pass $pass1$suff e $pass2$suff"
    cd $bam_dir
    samtools merge $all_name$suff $pass1$suff $pass2$suff

    echo "Sort per nome read del file $all_name$suff"
    samtools sort -o $all_sorted$name$suff -n -T $temp$suff $all_name$suff

    echo "Sort per posizione read e creazione indice del file $all_name$suff"
    samtools sort -o $all_sorted$suff -T $temp$suff $all_name$suff
    samtools index -b $all_sorted$suff
    cd ..
}

# ordina i file pass creati e crea l'indice
sort_bam() {
    suff=".bam"
    unique="unique_reads"
    single="single_reads"
    multiple="multiple_reads"
    SUFF="_sorted"
    temp="temp"

    cd $ris_dir
    for read in $unique $single $multiple; do
        samtools sort -o $read$SUFF$suff -T $temp$suff $read$suff
        samtools index -b $read$SUFF$suff
    done
    cd ..
}

# crea dai file wig i relativi file tdf
wig_to_tdf() {
    echo "Creazione file tdf"
    suff="wig"
    SUFF="tdf"
    genome="reference_A_laidlawii/reference_A_laidlawii.fasta.fai"

    cd $ris_dir
    for s in $(ls *.$suff); do
        igvtools toTDF $s $s.$SUFF $genome
    done
    cd ..
}


# genera i grafici con gnuplot
plot() {
    echo "Generazione grafici lunghezza inserti"
    graph="plot.gnuplot"
    gnuplot $graph
}

# main
while getopts ":a:b:gms:tw" flag; do
    case $flag in
        a )
            merge_pass_sort
            start_reseq $OPTARG
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
        s ) start_reseq $OPTARG
            ;;
        t ) sort_bam
            ;;
        w ) wig_to_tdf
            ;;
        * ) 
            if [ $flag == ":" ]; then
                echo "Per eseguire la conversione dei file bam in sam è necessario inserire la cartella dove sono salvati i file."
            else
                echo "Opzione non presente, nessuna operazione da svolgere."
            fi
            ;;
    esac
done