

convert ../figures/plots_network_stats_two_types2.pdf -fill black \
-font Helvetica-bold -pointsize 20 \
label:' A ' \
-gravity northwest -geometry +4+4 \
-composite ../figures/plots_network_stats_two_types2_annotated.png

convert ../figures/plots_network_stats_within_type2.pdf -fill black \
-font Helvetica-bold -pointsize 20 \
label:' B ' \
-gravity northwest -geometry +4+4 \
-composite ../figures/plots_network_stats_within_type2_annotated.png


montage -mode concatenate -density 300 ../figures/plots_network_stats_two_types2_annotated.png ../figures/plots_network_stats_within_type2_annotated.png  ../figures/plots_network_stats_montage.png