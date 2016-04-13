## combine the network plots so that it is easier to put in the manuscript


 convert ../figures/network_LSA_gp23_to_gp23.graphml_viral_groups_coloured.png -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' A ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_gp23_to_gp23.graphml_viral_groups_coloured_annotated.png

## flipped and rotated so that the networks are in the same orientation


 convert ../figures/network_LSA_gp23_to_gp23.graphml_with_modules.png -flip -rotate 150 -trim  -fill black \
-font Helvetica-bold -pointsize 70 \
label:' B ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_gp23_to_gp23.graphml_with_modules_annotated.png


montage -mode concatenate -density 300 ../figures/network_LSA_gp23_to_gp23.graphml_viral_groups_coloured_annotated.png  ../figures/network_LSA_gp23_to_gp23.graphml_with_modules_annotated.png  ../figures/network_LSA_gp23_to_gp23.graphml_montage.png


convert ../figures/network_LSA_MPL_to_MPL.graphml_viral_groups_coloured.png -flip -rotate 60 -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' C ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_MPL_to_MPL.graphml_viral_groups_coloured_annotated.png


convert ../figures/network_LSA_MPL_to_MPL.graphml_with_modules.png -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' D ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_MPL_to_MPL.graphml_with_modules_annotated.png


montage -tile 2x2 -geometry +5+5 -density 300 ../figures/network_LSA_gp23_to_gp23.graphml_viral_groups_coloured_annotated.png  ../figures/network_LSA_gp23_to_gp23.graphml_with_modules_annotated.png ../figures/network_LSA_MPL_to_MPL.graphml_viral_groups_coloured_annotated.png ../figures/network_LSA_MPL_to_MPL.graphml_with_modules_annotated.png ../figures/network_LSA_gp23_to_gp23.graphml_montage.png
#montage -mode concatenate -density 300 ../figures/network_LSA_MPL_to_MPL.graphml_viral_groups_coloured_annotated.png ../figures/network_LSA_MPL_to_MPL.graphml_with_modules_annotated.png ../figures/network_LSA_MPL_to_MPL.graphml_montage.png



convert ../figures/network_LSA_euk_to_bac.graphml.png -fill black \
-font Helvetica-bold -pointsize 70 \
label:' A ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_euk_to_bac.graphml_annotated.png

convert ../figures/network_LSA_euk_to_bac.graphml_with_modules.png -rotate 220 -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' B ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_euk_to_bac.graphml_with_modules_annotated.png


montage -mode concatenate -density 300 ../figures/network_LSA_euk_to_bac.graphml_annotated.png ../figures/network_LSA_euk_to_bac.graphml_with_modules_annotated.png  ../figures/network_LSA_euk_to_bac.graphml_montage.png



convert ../figures/network_LSA_bac_to_gp23.graphml_viral_groups_coloured.png -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' C ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_bac_to_gp23.graphml_viral_groups_coloured_annotated.png

convert ../figures/network_LSA_bac_to_gp23.graphml_with_modules.png  -rotate 110 -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' D ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_bac_to_gp23.graphml_with_modules_annotated.png


#montage -mode concatenate -density 300 ../figures/network_LSA_bac_to_gp23.graphml_viral_groups_coloured_annotated.png ../figures/network_LSA_bac_to_gp23.graphml_with_modules_annotated.png  ../figures/network_LSA_bac_to_gp23.graphml_montage.png



convert ../figures/network_LSA_euk_to_MPL.graphml_viral_groups_coloured.png -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' E ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_euk_to_MPL.graphml_viral_groups_coloured_annotated.png

convert ../figures/network_LSA_euk_to_MPL.graphml_with_modules.png -rotate 270 -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' F ' \
-gravity northwest -geometry +20+20 \
-composite ../figures/network_LSA_euk_to_MPL.graphml_with_modules_annotated.png


#montage -mode concatenate -density 300 ../figures/network_LSA_euk_to_MPL.graphml_viral_groups_coloured_annotated.png ../figures/network_LSA_euk_to_MPL.graphml_with_modules_annotated.png  ../figures/network_LSA_euk_to_MPL.graphml_montage.png


montage -mode concatenate -density 300 -geometry 1200x1200+5+5 -tile 2x3 ../figures/network_LSA_euk_to_bac.graphml_annotated.png ../figures/network_LSA_euk_to_bac.graphml_with_modules_annotated.png ../figures/network_LSA_bac_to_gp23.graphml_viral_groups_coloured_annotated.png ../figures/network_LSA_bac_to_gp23.graphml_with_modules_annotated.png ../figures/network_LSA_euk_to_MPL.graphml_viral_groups_coloured_annotated.png ../figures/network_LSA_euk_to_MPL.graphml_with_modules_annotated.png ../figures/network_LSA_euk_to_bac.graphml_montage.png


## do for bac and euks too

convert ../figures/network_LSA_bac_to_bac.graphml.png -fill black \
-font Helvetica-bold -pointsize 70 \
label:' A ' \
-gravity northwest -geometry +20+20 \
-composite -trim ../figures/network_LSA_bac_to_bac.graphml_annotated.png

convert ../figures/network_LSA_bac_to_bac.graphml_with_modules.png -flip -fill black \
-font Helvetica-bold -pointsize 70 \
label:' B ' \
-gravity northwest -geometry +20+20 \
-composite -trim ../figures/network_LSA_bac_to_bac.graphml_with_modules_annotated.png


#montage -mode concatenate -density 300 ../figures/network_LSA_bac_to_bac.graphml_annotated.png ../figures/network_LSA_bac_to_bac.graphml_with_modules_annotated.png  ../figures/network_LSA_bac_to_bac.graphml_montage.png




convert ../figures/network_LSA_euk_to_euk.graphml.png -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' C ' \
-gravity northwest -geometry +20+20 \
-composite -trim ../figures/network_LSA_euk_to_euk.graphml_annotated.png

convert ../figures/network_LSA_euk_to_euk.graphml_with_modules.png -flip -rotate 310 -trim -fill black \
-font Helvetica-bold -pointsize 70 \
label:' D ' \
-gravity northwest -geometry +20+20 \
-composite -trim ../figures/network_LSA_euk_to_euk.graphml_with_modules_annotated.png


montage -mode concatenate -density 300 ../figures/network_LSA_bac_to_bac.graphml_annotated.png ../figures/network_LSA_bac_to_bac.graphml_with_modules_annotated.png ../figures/network_LSA_euk_to_euk.graphml_annotated.png ../figures/network_LSA_euk_to_euk.graphml_with_modules_annotated.png   ../figures/network_LSA_euk_to_euk.graphml_montage.png
