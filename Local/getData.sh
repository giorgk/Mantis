#!/bin/bash
# This is working only for the xml files.
# For the tifs themself it seems that it downloads data but the data cannot be read by matlab.

# rm -r Ngw_*
for ((year=1945; year<=2050;year=year+15)); do
	#curl -o Ngw_$year.tif https://github.com/ucd-cws/nitrates-cv/blob/master/$year/Ngw.tif
	curl -o Ngw_$year.tif.xml https://raw.githubusercontent.com/ucd-cws/nitrates-cv/master/1945/Ngw.tif.xml
done

#curl -o Ngw_1945.tif https://github.com/ucd-cws/nitrates-cv/blob/master/1945/Ngw.tif?raw=true
