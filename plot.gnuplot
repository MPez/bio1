set title "Genomic insert length"
set terminal pdf
set xlabel "Numero inserti"
set ylabel "Lunghezza inserti (bp)"
set pointsize 0.1
set xtics rotate

set output "risultati/insert_length1.pdf"
plot "risultati/insert_length.txt" every ::1::50000 using 1:3 with dots title "Inserti da 1 a 50000"

set output "risultati/insert_length1-1.pdf"
plot "risultati/insert_length.txt" every ::50001::100000 using 1:3 with dots title "Inserti da 50001 a 100000"

set output "risultati/insert_length2.pdf"
plot "risultati/insert_length.txt" every ::100001::200000 using 1:3 with dots title "Inserti da 100001 a 200000"

set output "risultati/insert_length3.pdf"
plot "risultati/insert_length.txt" every ::200001::250000 using 1:3 with dots title "Inserti da 200001 a 250000"

set output "risultati/insert_length3-1.pdf"
plot "risultati/insert_length.txt" every ::250001::300000 using 1:3 with dots title "Inserti da 250001 a 300000"

set output "risultati/insert_length4.pdf"
plot "risultati/insert_length.txt" every ::300001::350000 using 1:3 with dots title "Inserti da 300001 a 350000"

set output "risultati/insert_length4-1.pdf"
plot "risultati/insert_length.txt" every ::350001::400000 using 1:3 with dots title "Inserti da 350001 a 350000"

set output "risultati/insert_length5.pdf"
plot "risultati/insert_length.txt" every ::400001::500000 using 1:3 with dots title "Inserti da 400001 a 500000"

set output "risultati/insert_length6.pdf"
plot "risultati/insert_length.txt" every ::500001::550000 using 1:3 with dots title "Inserti da 500001 a 550000"

set output "risultati/insert_length6-1.pdf"
plot "risultati/insert_length.txt" every ::550001::600000 using 1:3 with dots title "Inserti da 550001 a 600000"

set output "risultati/insert_length7.pdf"
plot "risultati/insert_length.txt" every ::600001 using 1:3 with dots title "Inserti da 600001"

set output "risultati/discarded_insert_length.pdf"
plot "risultati/discarded_insert_length.txt" using 1:3 with dots title "Inserti scartati"