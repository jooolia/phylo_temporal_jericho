
## convert pdfs to pngs for use in presentations

for file in ../figures/*.pdf
	do
	echo $file
	short_name=`basename -s ".pdf" $file`
	convert -density 300 -trim $file ../figures/$short_name".png"
	done