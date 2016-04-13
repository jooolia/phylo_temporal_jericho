

## pull down chlorophyll files and run the rscript on them to generate plots. 


June20=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011171.L3m_DAY_CHL_chlor_a_4km.bz2
June21=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011172.L3m_DAY_CHL_chlor_a_4km.bz2
June22=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011173.L3m_DAY_CHL_chlor_a_4km.bz2
June23=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011174.L3m_DAY_CHL_chlor_a_4km.bz2
June24=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011175.L3m_DAY_CHL_chlor_a_4km.bz2
June25=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011176.L3m_DAY_CHL_chlor_a_4km.bz2
June26=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011177.L3m_DAY_CHL_chlor_a_4km.bz2
June27=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011178.L3m_DAY_CHL_chlor_a_4km.bz2
June28=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011179.L3m_DAY_CHL_chlor_a_4km.bz2
June29=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011180.L3m_DAY_CHL_chlor_a_4km.bz2
June30=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011181.L3m_DAY_CHL_chlor_a_4km.bz2
July01=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011182.L3m_DAY_CHL_chlor_a_4km.bz2
July02=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011183.L3m_DAY_CHL_chlor_a_4km.bz2
July03=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011184.L3m_DAY_CHL_chlor_a_4km.bz2
July04=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011185.L3m_DAY_CHL_chlor_a_4km.bz2
July05=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011186.L3m_DAY_CHL_chlor_a_4km.bz2
July06=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011187.L3m_DAY_CHL_chlor_a_4km.bz2
July07=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2011188.L3m_DAY_CHL_chlor_a_4km.bz2
June18_June25=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20111692011176.L3m_8D_CHL_chlor_a_4km.bz2
June26_July03=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20111772011184.L3m_8D_CHL_chlor_a_4km.bz2
July04_July11=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20111852011192.L3m_8D_CHL_chlor_a_4km.bz2
June2011=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20111522011181.L3m_MO_CHL_chlor_a_4km.bz2
June2010=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20101522010181.L3m_MO_CHL_chlor_a_4km.bz2
June2012=http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20121532012182.L3m_MO_CHL_chlor_a_4km.bz2

for files in $June20 $June21 $June22 $June23 $June24 $June25 $June26 $June27 $June28 $June29 $June30 $July01 $July02 $July03 $July04 $July05 $July06 $July07 $June18_June25 $June26_July03 $July04_July11 $June2011 $June2010 $June2012
#chlorophyll_a_files  
do
	echo $files

	saved=`basename $files` 
	unzipped=`basename -s ".bz2" $files`
	ncfile=`basename -s ".bz2" $files | cut -d. -f1`
	wget $files 
	mv $saved ../data/$saved
	bunzip2  ../data/$saved
	/home/julia/h4tonccf/h4tonccf_nc4 ../data/$unzipped

	## gives truncated version of file
	echo ../data/$ncfile".nc"
	Rscript 30_testing_chlorophyl_satellite_image.R ../data/$ncfile".nc" ../figures/$unzipped"_chla_plot.pdf"

done


 