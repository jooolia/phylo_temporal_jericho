## blast viral OTUs from the good lsa. 

seqtk=/home/julia/seqtk/seqtk

good_lsa_results=../results/LSA_tables/good_lsa.csv

miseq_gp23_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_gp23_pear.assembled_annotated_quality_length_filt_otus_95.00.fasta


grep gp23 $good_lsa_results | cut -d , -f 2  > ../results/LSA_tables/good_lsa_gp23_col.txt
grep gp23 $good_lsa_results | cut -d , -f 3  >> ../results/LSA_tables/good_lsa_gp23_col.txt

grep gp23 ../results/LSA_tables/good_lsa_gp23_col.txt | sort |uniq > ../results/LSA_tables/good_lsa_gp23.txt

sed -i 's/"//g' ../results/LSA_tables/good_lsa_gp23.txt
sed -i 's/gp23//g' ../results/LSA_tables/good_lsa_gp23.txt

$seqtk subseq $miseq_gp23_amplicons ../results/LSA_tables/good_lsa_gp23.txt > ../results/LSA_tables/good_lsa_gp23.fasta


## want to blast and get list of the top blast hits back..
## biopython!!!


python blast_good_lsa_otus_on_NCBI.py ../results/LSA_tables/good_lsa_gp23.fasta ../results/LSA_tables/good_lsa_gp23_blast_hits.csv


### also do for mpl
miseq_MPL_amplicons=../../Miseq_Initial_Run_Processing/results/Divided_by_primers/All_libs_combined/Translated/Picked_OTUs/Total_MPL_pear.assembled_lines_annotated_quality_length_filt_otus_95.00.fasta


grep MPL $good_lsa_results | cut -d , -f 2  > ../results/LSA_tables/good_lsa_MPL_col.txt
grep MPL $good_lsa_results | cut -d , -f 3  >> ../results/LSA_tables/good_lsa_MPL_col.txt

grep MPL ../results/LSA_tables/good_lsa_MPL_col.txt | sort |uniq > ../results/LSA_tables/good_lsa_MPL.txt

sed -i 's/"//g' ../results/LSA_tables/good_lsa_MPL.txt
sed -i 's/MPL//g' ../results/LSA_tables/good_lsa_MPL.txt

$seqtk subseq $miseq_MPL_amplicons ../results/LSA_tables/good_lsa_MPL.txt > ../results/LSA_tables/good_lsa_MPL.fasta


python blast_good_lsa_otus_on_NCBI.py ../results/LSA_tables/good_lsa_MPL.fasta ../results/LSA_tables/good_lsa_MPL_blast_hits.csv

