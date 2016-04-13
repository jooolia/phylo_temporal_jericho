

## Raxml


## with bootstrapping -D is seach convergence criterion. Can turn off for final trees?

­-f c check if the alignment can be properly read by RAxML 
raxmlHPC ­m PROTGAMMAAUTO  ## model test...
 ­-o Rat outgroup.. multiple ­o Rat,Mouse,
­­--silent ## if do initial check with -f c this option could save time
­­--no­seq­check ## same as above

-T threads
raxmlHPC ­-b 12345 ­-p 333 ­-N 100 ­-s alg ­-m PROTGAMMAVT ­-n TEST ­



Now, if you want to run a full analysis, i.e., BS and ML search type:
raxmlHPC ­-f a ­-x 12345 ­-p 333 -N 100 ­-m PROTGAMMAVT ­-s ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta ­-n RdRptree 

/home/labop/RAxML-7.2.8-ALPHA/raxmlHPC-PTHREADS -T 8 -f -a ­-x 12345 ­-p 333 -N 100 ­-m PROTGAMMAVT ­-s ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta ­-n RdRptree

/Data/software/standard-RAxML-master/raxmlHPC-PTHREADS -f c -x 12345 -p 333 -N 100 -m PROTGAMMAVT -s ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta -n RdRptree


/Data/software/standard-RAxML-master/raxmlHPC-PTHREADS -T 8 -f a -x 12345 -N 100 -p 333 -m PROTGAMMAVT -s ../results/RdRp_95_miseq_data_with_env_iso_and_ref_align_clustal_edited.fasta -n RdRptree
