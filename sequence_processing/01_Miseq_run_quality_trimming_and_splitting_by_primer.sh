## Overall Miseq run processing from beginning to the end analysis.
## Author: Julia Gustavsen
## Last updated: 10 December 2014

# To-do: 
# what was the qiime script used to separate the libraries?
# Add back in that I want to caluclate the length of the sequences, after trimming and when they are divided by primers.
## need to run this before nounset or else get error. 
source /home/labop/Data/Julia/Julia_python/bin/activate
## Purpose: Initial processing of sequences. 
#
# Quality checking of libraries, trim off adapters, divide by primers, concatenate libraries together
## This was too stringent for the script. 
set -o nounset -o pipefail -o errexit

## set data directory so it is easy to move around
data_dir=../data


### Set up python virtual env on server scr
## should probably do the same on my laptop. 

echo "############ Unzip the data files ##################"

 # for files in $data_dir/*.fastq.gz
 #   do
 #       echo $files
 #       gunzip $files
 #    done



############### Qiime for libs that need to be split by barcode ################

## run the file to split the libraries after running the qiime script
# bash grep_to_divide_by_lib_from_qiime_split_lib.sh

################ transfer files via filezilla to my laptop.


echo "##################### Rename Libs #########################"

### Will add in the mapping renaming python file here.
 #cp ../../Data_Illumina_Original/*.fastq ../data/


## uses files in the ../data directory
## Rename files using mapping file which has all of the VC information in it.

echo "Rename libraries to something based on my metadata file"
 python renaming_library_files_according_to_mapping_file.py

echo "################ Fast QC on the initial unadultered files ##############"
 Fastq_stats_dir=../results/Fastq_stats
 Fastq_qual_stats_dir=../results/Fastq_stats/Fastq_qual
 fastqc_dir=../results/FastQC
 Primer_counts=../results/Primer_counts
 mkdir -p $Fastq_stats_dir
 mkdir -p $Fastq_qual_stats_dir
 mkdir -p $fastqc_dir
 mkdir -p $Primer_counts

  for files in $data_dir/*.fastq
  do
    echo $files	
 #   /home/julia/FastQC/fastqc $files -o $fastqc_dir
    # Fixing the header to undo what qiime did to it.
    sed -i '/^@/s/ orig_bc=.*bc_diffs=0//g;/^@/s/ /_/g' $files 
    echo "Count Primers in the files"
    ## usage: python Counting_primers_with_python.py file.fastq
    python Counting_primers_with_python.py $files $Primer_counts
  done

#rm $fastqc_dir/*_fastqc.zip



results_dir=../results

echo "###### Trimmomatic ######"   
 ## Not doing quality trimming before merging 
 ## Trim sequencing adapters from Illumina reads
 Adapters=../data/Adapters_Julia_Sep_2014.fa

 file=../../JerichoAndSOGsequencing/data/Primers_reg_ex_Julia_July2013.txt
 Primer_divided_dir=../results/Divided_by_primers/
 Trimmomatic=/home/labop/Trimmomatic-0.32/trimmomatic-0.32.jar
 fastx_quality=/home/labop/fastx_bin/fastx_quality_stats
 Length_dir=../results/Length_data

  mkdir -p $Primer_divided_dir
  mkdir -p $Length_dir


  for number in {1..16} {18..59} 98 99 16redo 17redo

  do
   FastQ_R1_orig=$data_dir/"Lib_pool"$number"_R1.fastq"
   FastQ_R2_orig=$data_dir/"Lib_pool"$number"_R2.fastq"


   FastQ_R1_orig_quality_stats=$Fastq_qual_stats_dir/"Fastq_quality_stats_Lib_pool"$number"_R1.txt"
   FastQ_R2_orig_quality_stats=$Fastq_qual_stats_dir/"Fastq_quality_stats_Lib_pool"$number"_R2.txt"

   FastQ_R1_out_paired=../results/"Lib_pool"$number"_R1_trimmed_paired.fastq"
   FastQ_R1_out_unpaired=../results/"Lib_pool"$number"_R1_trimmed_unpaired.fastq"
   FastQ_R2_out_paired=../results/"Lib_pool"$number"_R2_trimmed_paired.fastq"
   FastQ_R2_out_unpaired=../results/"Lib_pool"$number"_R2_trimmed_unpaired.fastq"


   FastQ_R1_fastq_quality_stats=$Fastq_qual_stats_dir/"Fastq_quality_stats_Lib_pool"$number"_R1_trimmed_paired.txt"
   FastQ_R2_fastq_quality_stats=$Fastq_qual_stats_dir/"Fastq_quality_stats_Lib_pool"$number"_R2_trimmed_paired.txt"
   ### Trim the sequences using the adapters####
     # Keeping minlength as 150 so that I still get high quality 18s sequences. 
     #  Just want to trim the adapters, but not get rid of adapters so that I can use uparse pipeline. 
    java -jar $Trimmomatic PE -phred33 -threads 8 -trimlog trimmomatic_test.log $FastQ_R1_orig $FastQ_R2_orig $FastQ_R1_out_paired $FastQ_R1_out_unpaired $FastQ_R2_out_paired $FastQ_R2_out_unpaired ILLUMINACLIP:$Adapters:2:30:10 MINLEN:150

  #  /home/julia/FastQC/fastqc $FastQ_R1_out_paired -o $fastqc_dir
  #  /home/julia/FastQC/fastqc $FastQ_R2_out_paired -o $fastqc_dir

  # use robert edgar fastq stat on untrimmed data  


$fastx_quality -Q33 -i $FastQ_R1_orig -o $FastQ_R1_orig_quality_stats
$fastx_quality -Q33 -i $FastQ_R2_orig -o $FastQ_R2_orig_quality_stats
$fastx_quality -Q33 -i $FastQ_R1_out_paired -o $FastQ_R1_fastq_quality_stats
$fastx_quality -Q33 -i $FastQ_R2_out_paired -o $FastQ_R2_fastq_quality_stats



# Count the length of the sequences
   python count_sequence_length.py $FastQ_R1_out_paired $Length_dir
   python count_sequence_length.py $FastQ_R2_out_paired $Length_dir

   #   check primer counts

   echo "########## Split libraries by primer ######################"


  python split_trimmed_libs_by_primer.py $FastQ_R1_out_paired $FastQ_R2_out_paired

    done



##Remove the fastqc zip files to clean-up
 #rm $fastqc_dir/*_fastqc.zip

echo "################### Zip up the Fastq Files #################"

#gzip $data_dir/*.fastq


echo "######################### Merge reads by primer ##############"
## Merging is better with Pear if I do trimming of adapters but not quality before hand.
## Want to merge files that look like Lib_pool98_16s_R1.fastq with Lib_pool98_16s_R2.fastq

pear_program=/home/labop/PEAR/bin/pear-0.9.6-bin-64
# old /home/labop/pear-0.9.5-bin-64/pear-0.9.5-64



mkdir -p ../results/Merge_counts/

echo "****** Merging reads together by primer *******"

## for testing let some errors go by...the files might not be there and thus should not fail. 
 set +u 
 set +e
for number in {1..16} {18..59} 98 99 16redo 17redo
do
  for primer in MPL 18s gp23 
  do 
    Primer_lib_R1=$Primer_divided_dir"Lib_pool"$number"_"$primer"_R1.fastq"
    Primer_lib_R2=$Primer_divided_dir"Lib_pool"$number"_"$primer"_R2.fastq"
   Merged_Fastq=$Primer_divided_dir"Lib_pool"$number"_"$primer"_pear"
    echo $Primer_lib_R1
    echo $Primer_lib_R2
    echo $Merged_Fastq

       if [ "$primer" = "MPL" ]
    	then
    	min_length=400
    	max_length=550

    # Alvin is using these settings: -n 251 -u 0 -q 20 -t 5 -j 12
    # I could change the settings for each primer! 400 for MPL, 16s, gp23, 200 for 18s
    # -u 0 will give more stringent settings because it is proportion of uncalled reads
    #-j 8 threads on this labpto
    # use setting w different than 1 for the short overlaps(like the 16s??)
    # could set min -n and max length -m for sequences...maybe check histograms first
    # -g 2 changes the statistical test to recover more reads with short overlaps. Could try for 16s?
    # I think I want to get rid of reads with "N"s in them


    $pear_program -f $Primer_lib_R1 -r $Primer_lib_R2 -o $Merged_Fastq -v 20 -n $min_length -m $max_length -j 8 -g 2 

elif [ "$primer" = "18s" ]
 then
    	min_length=150
    	max_length=250
    	
   $pear_program -f $Primer_lib_R1 -r $Primer_lib_R2 -o $Merged_Fastq -v 20 -n $min_length -m $max_length -u 0 -j 8 

elif [ "$primer" = "gp23" ]
 then
    	min_length=400
    	max_length=600
    	
   $pear_program -f $Primer_lib_R1 -r $Primer_lib_R2 -o $Merged_Fastq -v 20 -n $min_length -m $max_length -u 0 -j 8 -g 2 
fi
 Merged_Fastq_assembled=$Primer_divided_dir"Lib_pool"$number"_"$primer"_pear.assembled.fastq"
    python count_merged_primer_files.py $Primer_lib_R1 $Primer_lib_R2 $Merged_Fastq_assembled
   done
done     


 set -u 
 set -e
# clean up a bit
 rm ../results/Divided_by_primers/*unassembled*


echo "############### Concatenate sequences from all libraries together by primer ####################"


## To hold quality stats


Fastq_qual_by_primer=../results/Fastq_stats/Fastq_qual_by_primer/
Length_dir=../results/Length_data
mkdir -p $Fastq_qual_by_primer
mkdir -p $Length_dir


## temporary
set +u +e
## First need to add the library identifyers to all the sequences
for number in {1..16} {18..59} 98 99 16redo 17redo
do
  for primer in MPL 18s gp23 
  do 

  Merged_Fastq_=$Primer_divided_dir"Lib_pool"$number"_"$primer"_pear.assembled.fastq"
  Merged_Fastx_qual_=$Fastq_qual_by_primer"Lib_pool"$number"_"$primer"_pear.assembled_fastx_qual.txt"
  Merged_Fastq_lines_annotated=$Primer_divided_dir"Lib_pool"$number"_"$primer"_pear.assembled_lines_annotated.fastq"
  Merged_Fastq_lines_annotated_quality_trim=$Primer_divided_dir"Lib_pool"$number"_"$primer"_pear.assembled_lines_annotated_quality_trim.fastq" 
  Merged_Fastx_qual_trimmed=$Fastq_qual_by_primer"Lib_pool"$number"_"$primer"_pear.assembled_quality_trimmed_fastx_qual.txt"



      python count_sequence_length.py $Merged_Fastq_ $Length_dir

      python annotate_lines_of_fastq_libs.py $Merged_Fastq_ $Merged_Fastq_lines_annotated 


      $fastx_quality -Q33 -i $Merged_Fastq_ -o $Merged_Fastx_qual_

     java -jar $Trimmomatic SE -threads 8 -phred33  $Merged_Fastq_lines_annotated $Merged_Fastq_lines_annotated_quality_trim LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:175

      $fastx_quality -Q33 -i $Merged_Fastq_lines_annotated_quality_trim -o $Merged_Fastx_qual_trimmed

  done
done
set -u -e
## need to combine everything for each primer set
mkdir -p ../results/Divided_by_primers/All_libs_combined


  for primer in MPL 18s gp23 
  do
    cat ../results/Divided_by_primers/"Lib_pool"*"_"$primer"_pear.assembled_lines_annotated_quality_trim.fastq" > ../results/Divided_by_primers/All_libs_combined/"Total_"$primer"_pear.assembled_lines_annotated_quality_trim.fastq"
done



#rm ../results/Divided_by_primers/"Lib_pool"*"_"*"_pear.discarded.fastq"
#rm  ../results/Divided_by_primers/"Lib_pool"*"_"*"_pear.assembled_lines_annotated_quality_trim.fastq"




set +u +e
echo "#################### 16S that were not expected to merge and even those that were expected to merge..#####################"
## for 16s
## for libraries that I did not attempt to merge
for number in {1..16} {18..59} 98 99 16redo 17redo
do

Fastq_R1_=$Primer_divided_dir"Lib_pool"$number"_16s_R1.fastq"
Fastq_R2_=$Primer_divided_dir"Lib_pool"$number"_16s_R2.fastq"

Fastx_qual_Fastq_R1_=$Fastq_qual_by_primer"Lib_pool"$number"_16s_R1_fastx_qual.txt"
Fastx_qual_Fastq_R2_=$Fastq_qual_by_primer"Lib_pool"$number"_16s_R2_fastx_qual.txt"


Fastq_R1_lines_annotated=$Primer_divided_dir"Lib_pool"$number"_16s_R1_annotated.fastq"
Fastq_R2_lines_annotated=$Primer_divided_dir"Lib_pool"$number"_16s_R2_annotated.fastq"


Fastq_R1_lines_annotated_trimmed_paired=$Primer_divided_dir"Lib_pool"$number"_16s_R1_annotated_trimmed_paired.fastq"
Fastq_R1_lines_annotated_trimmed_unpaired=$Primer_divided_dir"Lib_pool"$number"_16s_R1_annotated_trimmed_unpaired.fastq"
Fastq_R2_lines_annotated_trimmed_paired=$Primer_divided_dir"Lib_pool"$number"_16s_R2_annotated_trimmed_paired.fastq"
Fastq_R2_lines_annotated_trimmed_unpaired=$Primer_divided_dir"Lib_pool"$number"_16s_R2_annotated_trimmed_unpaired.fastq"


Fastx_qual_Fastq_R1_trimmed_paired=$Fastq_qual_by_primer"Lib_pool"$number"_16s_R1_trimmed_paired_fastx_qual.txt"
Fastx_qual_Fastq_R2_trimmed_paired=$Fastq_qual_by_primer"Lib_pool"$number"_16s_R2_trimmed_paired_fastx_qual.txt"

python count_sequence_length.py $Fastq_R1_ $Length_dir
python count_sequence_length.py $Fastq_R2_ $Length_dir
 

 $fastx_quality -Q33 -i $Fastq_R1_ -o $Fastx_qual_Fastq_R1_
 $fastx_quality -Q33 -i $Fastq_R2_ -o $Fastx_qual_Fastq_R2_
 
python annotate_lines_of_fastq_libs.py $Fastq_R1_ $Fastq_R1_lines_annotated 
python annotate_lines_of_fastq_libs.py $Fastq_R2_ $Fastq_R2_lines_annotated

java -jar $Trimmomatic PE -threads 8 -phred33 -trimlog trimmomatic_test.log $Fastq_R1_lines_annotated $Fastq_R2_lines_annotated $Fastq_R1_lines_annotated_trimmed_paired $Fastq_R1_lines_annotated_trimmed_unpaired $Fastq_R2_lines_annotated_trimmed_paired $Fastq_R2_lines_annotated_trimmed_unpaired LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:175

$fastx_quality -Q33 -i $Fastq_R1_lines_annotated_trimmed_paired -o $Fastx_qual_Fastq_R1_trimmed_paired
$fastx_quality -Q33 -i $Fastq_R2_lines_annotated_trimmed_paired -o $Fastx_qual_Fastq_R2_trimmed_paired


done
set -u -e

for read in 1 2 
do
    cat ../results/Divided_by_primers/"Lib_pool"*"_16s_R"$read"_annotated_trimmed_paired.fastq" > ../results/Divided_by_primers/All_libs_combined/"Total_16s_R"$read"_annotated_trimmed_paired.fastq"
done

## clean up
#rm ../results/Divided_by_primers/"Lib_pool"*"_16s_R"*"_annotated_trimmed_paired.fastq"
#rm $Primer_divided_dir"Lib_pool"*"annotated_trimmed_unpaired.fastq"

#### So for 16s I end up with 4 files with all the libs:
## Total_16s_R1.fastq
## Total_16s_R2.fastq





####### Look at overall quality of the libraries ##############


### This seems out of place here...
## Want to copy all the FastQC summary files by different type to new dir so I can process easily in R.

# mkdir -p ../results/FastQC/Collected_from_all_libs

# #This copies the files to a new directory and renames them after their directory. Cool
# for files in ../results/FastQC/*/fastqc_data.txt
# do
#  parent_dir_name=`echo $files | cut -d \/ -f4`

# cp $files ../results/FastQC/Collected_from_all_libs/$parent_dir_name"_data.txt"
# done


# get rid of these empty files before reading into R. 
find $Fastq_qual_stats_dir -size 0 -delete
find $Fastq_qual_by_primer -size 0 -delete


### Now process in R! 
################## run files FastQC_plotting_all_libs_together.R and the markdown file

### clean up 

## do that soon!
#gzip -f ../results/*.fastq 
#gzip -f $Primer_divided_dir/*.fastq

set +u +e
deactivate
set -u -e