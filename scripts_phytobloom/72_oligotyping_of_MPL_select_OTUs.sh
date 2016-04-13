## entropy analysis of MPL


## get top 10 OTUs, pull out reads, derep, align and then edit those alignments


source /home/labop/Data/Julia/Julia_python/bin/activate


set -o nounset -o pipefail -o errexit


## Programs
## Programs

#usearch8=/home/julia/usearch/usearch8.0.1517_i86linux32
usearch8=/home/labop/usearch8/usearch8.0.1517_i86linux32

clustal_omega=/home/labop/clustal_omega/clustalo-1.2.0-Ubuntu-x86_64
#clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64
hmm_build=hmmbuild
#hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 

#seqtk=/home/julia/seqtk/seqtk
seqtk=seqtk





# #head $normalized_otu_table_top_10 -n 1 | sed s/\"//g > top_10_OTUs.txt
# 	for OTU in OTU_1 OTU_544 OTU_283 OTU_1588
# 	do

# 	oligotype ../results/MPL_"$OTU"like_reads_aligned_new_names.fasta" ../results/MPL"_"$OTU"like_reads_aligned_new_names.fasta-ENTROPY" -c 1 -A 10 -M 50 -o ../results/"oligotypingMPL"_"$OTU"like_reads"/ --gen-html 

# 		#echo "done with aligning dereplicated seqs for " $MPL
# 	done
# done


oligotype ../results/MPL_OTU_1like_reads_aligned_new_names.fasta ../results/MPL_OTU_1like_reads_aligned_new_names.fasta-ENTROPY -c 1 -A 10 -M 50 -o ../results/oligotypingMPL_OTU_1like_reads/ --gen-html 

oligotype ../results/MPL_OTU_544like_reads_aligned_new_names.fasta ../results/MPL_OTU_544like_reads_aligned_new_names.fasta-ENTROPY -c 1 -A 10 -M 50 -o ../results/oligotypingMPL_OTU_544like_reads/ --gen-html 

oligotype ../results/MPL_OTU_283like_reads_aligned_new_names.fasta ../results/MPL_OTU_283like_reads_aligned_new_names.fasta-ENTROPY -c 1 -A 10 -M 50 -o ../results/oligotypingMPL_OTU_283like_reads/ --gen-html

oligotype ../results/MPL_OTU_1588like_reads_aligned_new_names.fasta ../results/MPL_OTU_1588like_reads_aligned_new_names.fasta-ENTROPY -c 1 -A 10 -M 50 -o ../results/oligotypingMPL_OTU_1588like_reads/ --gen-html


oligotype ../results/MPL_OTU_1like_reads_aligned_all_times_new_names.fasta ../results/MPL_OTU_1like_reads_aligned_all_times_new_names.fasta-ENTROPY -c 1 -A 10 -M 50 -o ../results/oligotypingMPL_OTU_1like_reads_all_times/ --gen-html 

oligotype ../results/MPL_OTU_544like_reads_aligned_all_times_new_names.fasta ../results/MPL_OTU_544like_reads_aligned_all_times_new_names.fasta-ENTROPY -c 1 -A 10 -M 50 -o ../results/oligotypingMPL_OTU_544like_reads_all_times/ --gen-html 

oligotype ../results/MPL_OTU_283like_reads_aligned_all_times_new_names.fasta ../results/MPL_OTU_283like_reads_aligned_all_times_new_names.fasta-ENTROPY -c 1 -A 10 -M 50 -o ../results/oligotypingMPL_OTU_283like_reads_all_times/ --gen-html

oligotype ../results/MPL_OTU_1588like_reads_aligned_all_times_new_names.fasta ../results/MPL_OTU_1588like_reads_aligned_all_times_new_names.fasta-ENTROPY -c 1 -A 10 -M 50 -o ../results/oligotypingMPL_OTU_1588like_reads_all_times/ --gen-html


#python send_email_python.py $0

set +u +e
deactivate