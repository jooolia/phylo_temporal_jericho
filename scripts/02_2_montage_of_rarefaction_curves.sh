## montage of rarefaction curves



  convert ../figures/Rarefaction_curves_18S_original002.png -fill black \
-font Helvetica-bold -pointsize 28 \
label:' A ' \
-gravity northwest -geometry +10+12 \
-composite ../figures/Rarefaction_curves_18S_original_annotated.png

  convert ../figures/Rarefaction_curves_16S_R1_original002.png -fill black \
-font Helvetica-bold -pointsize 28 \
label:' B ' \
-gravity northwest -geometry +10+12 \
-composite ../figures/Rarefaction_curves_16S_R1_original_annotated.png

  convert ../figures/Rarefaction_curves_gp23_original002.png -fill black \
-font Helvetica-bold -pointsize 28 \
label:' C ' \
-gravity northwest -geometry +10+12 \
-composite ../figures/Rarefaction_curves_gp23_original_annotated.png

  convert ../figures/Rarefaction_curves_MPL_original002.png -fill black \
-font Helvetica-bold -pointsize 28 \
label:' D ' \
-gravity northwest -geometry +10+12 \
-composite ../figures/Rarefaction_curves_MPL_original_annotated.png




montage -mode concatenate -density 300 ../figures/Rarefaction_curves_18S_original_annotated.png ../figures/Rarefaction_curves_16S_R1_original_annotated.png  ../figures/Rarefaction_curves_gp23_original_annotated.png ../figures/Rarefaction_curves_MPL_original_annotated.png ../figures/Rarefaction_montage.png


