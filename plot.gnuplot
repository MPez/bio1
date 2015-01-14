set title "Genomic insert length"
set terminal pdf
set xlabel "Numero inserti"
set ylabel "Lunghezza inserti (bp)"
set pointsize 0.1

set output "insert_length1.pdf"
plot "insert_length.txt" every ::1::100000 using 1:3 with dots title "Inserti da 1 a 100000"

set output "insert_length2.pdf"
plot "insert_length.txt" every ::100001::200000 using 1:3 with dots title "Inserti da 100001 a 200000"

set output "insert_length3.pdf"
plot "insert_length.txt" every ::200001::300000 using 1:3 with dots title "Inserti da 200001 a 300000"

set output "insert_length4.pdf"
plot "insert_length.txt" every ::300001::400000 using 1:3 with dots title "Inserti da 300001 a 400000"

set output "insert_length5.pdf"
plot "insert_length.txt" every ::400001::500000 using 1:3 with dots title "Inserti da 400001 a 500000"

set output "insert_length6.pdf"
plot "insert_length.txt" every ::500001::600000 using 1:3 with dots title "Inserti da 500001 a 600000"

set output "insert_length7.pdf"
plot "insert_length.txt" every ::600001 using 1:3 with dots title "Inserti da 600001"

set output "discarded_insert_length.pdf"
plot "discarded_insert_length.txt" using 1:3 with dots title "Inserti scartati"