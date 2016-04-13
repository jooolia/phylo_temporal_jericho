### make the plots over time into little animationss

## test out with the 18s and MPL, AVS one right now

## names: cytoscape_network_AVS_MPL_and_18s-time_1.png



convert --adjoin `ls -1v ../figures/cytoscape_network_AVS_MPL_and_18s-time_*.png` ../figures/cytoscape_network_AVS_MPL_and_18s-time_plot.pdf



## the -lv puts them in the correct over
convert -delay 70 `ls -1v ../figures/cytoscape_network_AVS_MPL_and_18s-time_*.png` ../figures/cytoscape_network_AVS_MPL_and_18s-time.gif



convert -delay 70 `ls -1v ../figures/cytoscape_network_all-time_*.png` ../figures/cytoscape_network_all-time.gif