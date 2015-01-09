set title "Genomic insert length"
set terminal pdf
set xlabel "Numero inserti"
set ylabel "Lunghezza inserti (bp)"
set pointsize 0.1

set output "insert_length1.pdf"
plot "insert.txt" every ::1::300000 using 1:3 with dots title "Inserti da 1 a 300000"

set output "insert_length2.pdf"
plot "insert.txt" every ::300001 using 1:3 with dots title "Inserti da 300001"

set output "discarded_insert_length.pdf"
plot "discarded_insert.txt" using 1:3 with dots title "Inserti scartati"