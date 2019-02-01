## migth be useful to work on sequence alignment with all of the data? 
## FOR PROTEINS
## set up some directories for keeping downloaded references


## 
### Get references from Genbank
### Nucleotide alignments for viruses
source /home/labop/Data/Julia/Julia_python/bin/activate


set -o nounset -o pipefail -o errexit

## Programs
Sequence_directory=../results/Divided_by_primers/All_libs_combined
Translated_seqs_dir=$Sequence_directory/Translated

#clustal_omega=/home/julia/clustal/clustalo-1.2.0-Ubuntu-x86_64
clustal_omega=/home/labop/clustal_omega/clustalo-1.2.0-Ubuntu-x86_64

#hmm_build=/home/julia/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmbuild 
hmm_build=hmmbuild


## NCBI protein alignment

gp23_alignment=../results/Divided_by_primers/All_libs_combined/Translated/environmental_gp23_seqs_aligned_with_references_and_isolates_clustal_edited_realigned_edited.fasta
Sequences_gp23=$Translated_seqs_dir/Total_gp23_pear.assembled_lines_annotated_quality_trim.good.ng_good_by_blast.faa 


$hmm_build --informat afa gp23.hmm $gp23_alignment



echo "================================================"
echo "Split up file"
echo "================================================"			
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $Sequences_gp23 > $Sequences_gp23"_1L.fasta"

split -l 20000 -d -a 3 $Sequences_gp23"_1L.fasta" $Sequences_gp23"_1L.fasta"
# gives $Translated_seqs_dir/$basefilename"_clustal.fasta"[0-9][0-9]

for split_file in $Sequences_gp23"_1L.fasta"[0-9][0-9][0-9]
do

echo $split_file 


$clustal_omega -i $split_file  --hmm-in=gp23.hmm --outfmt=fasta -v --force -o $split_file"_clustal_aligned.fasta" 

done


## need to do profile alignment for 00 and up until there are no more...
 #$Sequences_gp23"_1L.fasta"[0-9][0-9]_clustal_aligned.fasta"

#cat $Translated_seqs_dir/$basefilename"derep_1L.fasta"[0-9][0-9]"_clustal_filtered_realigned_filtered_realigned_filtered.fasta" > $Translated_seqs_dir/$basefilename"_clustal_filtered_realigned_filtered_realigned_filtered.fasta" 



## need to do profile alignment for 00 and up until there are no more...

mv $Sequences_gp23"_1L.fasta000_clustal_aligned.fasta" $Sequences_gp23"initial_profile.fasta"
mv $Sequences_gp23"_1L.fasta001_clustal_aligned.fasta" $Sequences_gp23"initial_profile2.fasta"

$clustal_omega --p1=$Sequences_gp23"initial_profile.fasta"  --p2=$Sequences_gp23"initial_profile2.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profile.fasta"



## too slow and memory intensive...
for files in  $Sequences_gp23"_1L.fasta"{002..113}"_clustal_aligned.fasta"
do 
echo $files

$clustal_omega --p1=$files --p2=$Sequences_gp23"_clustal_aligned_all_profile.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profile.fasta"
done


## 114

mv $Sequences_gp23"_1L.fasta151_clustal_aligned.fasta" $Sequences_gp23"initial_profile3.fasta"
mv $Sequences_gp23"_1L.fasta152_clustal_aligned.fasta" $Sequences_gp23"initial_profile4.fasta"


$clustal_omega --p1=$Sequences_gp23"initial_profile3.fasta"  --p2=$Sequences_gp23"initial_profile4.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profileB.fasta"
for files in $Sequences_gp23"_1L.fasta"{153..200}"_clustal_aligned.fasta" 
do 
echo $files
$clustal_omega --p1=$files --p2=$Sequences_gp23"_clustal_aligned_all_profileB.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profileB.fasta"
done

### 201-250

mv $Sequences_gp23"_1L.fasta201_clustal_aligned.fasta" $Sequences_gp23"initial_profile5.fasta"
mv $Sequences_gp23"_1L.fasta202_clustal_aligned.fasta" $Sequences_gp23"initial_profile6.fasta"

$clustal_omega --p1=$Sequences_gp23"initial_profile5.fasta"  --p2=$Sequences_gp23"initial_profile6.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profileC.fasta"
for files in $Sequences_gp23"_1L.fasta"{203..250}"_clustal_aligned.fasta" 
do 
echo $files
$clustal_omega --p1=$files --p2=$Sequences_gp23"_clustal_aligned_all_profileC.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profileC.fasta"
done


## 250 + 

mv $Sequences_gp23"_1L.fasta251_clustal_aligned.fasta" $Sequences_gp23"initial_profile7.fasta"
mv $Sequences_gp23"_1L.fasta252_clustal_aligned.fasta" $Sequences_gp23"initial_profile8.fasta"


$clustal_omega --p1=$Sequences_gp23"initial_profile7.fasta"  --p2=$Sequences_gp23"initial_profile8.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profileD.fasta"
for files in $Sequences_gp23"_1L.fasta"{253..261}"_clustal_aligned.fasta" 
do 
echo $files
$clustal_omega --p1=$files --p2=$Sequences_gp23"_clustal_aligned_all_profileD.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profileD.fasta"
done



### 114 - 150


mv $Sequences_gp23"_1L.fasta114_clustal_aligned.fasta" $Sequences_gp23"initial_profile9.fasta"
mv $Sequences_gp23"_1L.fasta115_clustal_aligned.fasta" $Sequences_gp23"initial_profile10.fasta"


$clustal_omega --p1=$Sequences_gp23"initial_profile9.fasta"  --p2=$Sequences_gp23"initial_profile10.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profileE.fasta"
for files in $Sequences_gp23"_1L.fasta"{116..150}"_clustal_aligned.fasta" 
do 
echo $files
$clustal_omega --p1=$files --p2=$Sequences_gp23"_clustal_aligned_all_profileE.fasta" --outfmt=fasta -v --force -o $Sequences_gp23"_clustal_aligned_all_profileE.fasta"
done

