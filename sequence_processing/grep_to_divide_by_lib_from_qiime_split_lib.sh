## Divide libraries from fastq into each library based on the header.

for number in {1..10}
do 
echo $number
grep Julia.$number"_" -A 3 --no-group-separator /run/media/julia/OS/Users/Julia.Julia-PC/Documents/School/Thesis/ThesisProjects/Data_Illumina_Original/seqs_R1.fastq > ../data/"Julia_lib_"$number"_R1.fastq"
done

for number in {1..10}
do 
echo $number
grep Julia.$number"_" -A 3 --no-group-separator /run/media/julia/OS/Users/Julia.Julia-PC/Documents/School/Thesis/ThesisProjects/Data_Illumina_Original/seqs_R2.fastq > ../data/"Julia_lib_"$number"_R2.fastq"
done
