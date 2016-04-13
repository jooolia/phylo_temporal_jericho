## want to make networks in cytoscape using a bit of a script. 

cytoscape=/opt/Cytoscape_v3.2.1/cytoscape.sh 
/home/julia/Cytoscape_v3.2.1/cytoscape.sh

#/home/julia/Cytoscape_v3.2.0/cytoscape.sh -S cytoscape_script_16s_18s.txt -N ../figures/network_cooccurence_spearman_merged_16s_and_18s_all_types_proportion_above0_01.graphml

## add in a type for each here??
#/home/julia/Cytoscape_v3.2.0/cytoscape.sh -S cytoscape_script_16s_and_gp23.txt -N ../figures/network_cooccurence_spearman_merged_16s_and_gp23_all_types_proportion_above0_01.graphml

#/home/julia/Cytoscape_v3.2.0/cytoscape.sh -S cytoscape_script_18s_and_MPL_AVS.txt -N ../figures/network_cooccurence_spearman_merged_18s_and_AVS_MPL_all_types_proportion_above0_01.graphml




## time informed possibility
 #$cytoscape -S cytoscape_script_16s_and_gp23.txt -N ../figures/network_cooccurence_spearman_merged_16s_and_gp23_proportion_above0_01_temporal_data.graphml

 #$cytoscape -S cytoscape_script_18s_and_MPL_AVS.txt -N ../figures/network_cooccurence_spearman_merged_18s_and_AVS_MPL_proportion_above0_01_temporal_data.graphml

# $cytoscape -S cytoscape_script_16s_18s.txt -N ../figures/network_cooccurence_spearman_merged_16s_and_18s_proportion_above0_01_temporal_data.graphml

#$cytoscape -S cytoscape_script_all.txt -N ../figures/network_cooccurence_spearman_merged_all_proportion_above0_01_temporal_data.graphml

### networks with just viruses

#$cytoscape -S cytoscape_script_18s.txt -N ../figures/network_cooccurence_spearman_18s_proportion_above0_01_temporal_data.graphml

#$cytoscape -S cytoscape_script_16s.txt -N ../figures/network_cooccurence_spearman_16s_proportion_above0_01_temporal_data.graphml

# $cytoscape -S cytoscape_script_AVS.txt -N ../figures/network_cooccurence_spearman_AVS_proportion_above0_01_temporal_data.graphml

# $cytoscape -S cytoscape_script_gp23.txt -N ../figures/network_cooccurence_spearman_gp23_proportion_above0_01_temporal_data.graphml

# $cytoscape -S cytoscape_script_MPL.txt -N ../figures/network_cooccurence_spearman_MPL_proportion_above0_01_temporal_data.graphml



## simplified

$cytoscape -S cytoscape_script_16s_and_gp23_simplified.txt -N ../figures/network_cooccurence_spearman_merged_16s_and_gp23__temporal_data_only_edges_btwn_unlike.graphml

$cytoscape -S cytoscape_script_18s_and_MPL_AVS_simplified.txt -N ../figures/network_cooccurence_spearman_merged_18s_and_AVS_MPL__temporal_data_only_edges_btwn_unlike.graphml

#$cytoscape -S cytoscape_script_16s_18s_simplified.txt -N ../figures/network_cooccurence_spearman_merged_16s_and_18s__temporal_data_only_edges_btwn_unlike.graphml

# $cytoscape -S cytoscape_script_all_simplified.txt -N ../figures/network_cooccurence_spearman_merged_all__temporal_data_only_edges_btwn_unlike.graphml


# for files in ../figures/cytoscape_network*.pdf
# do
# 	base=`basename -s .pdf $files`
# 	echo $base
# 	convert -density 300 $files $base".png" 
# done