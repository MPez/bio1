#!/bin/bash
#preelaborazione dei file BAM forniti

samtools sort -o pass_reads1_sorted_name.bam -n -T prefix.bam pass_reads1.bam
samtools sort -o pass_reads2_sorted_name.bam -n -T prefix.bam pass_reads2.bam

samtools index -b pass_reads1_sorted_name.bam
samtools index -b pass_reads2_sorted_name.bam