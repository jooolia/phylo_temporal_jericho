seqtk=/home/labop/Data/software/seqtk/seqtk
#seqtk=/home/julia/seqtk/seqtk

for files in ../data/*.fastq
    do
    echo $files
    basename_reads=$(basename $files .fastq)
   $seqtk sample $files 10000 > ../data_small_test2/$basename_reads".fastq"
 echo "done subsampling"
   done


