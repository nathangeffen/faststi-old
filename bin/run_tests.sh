#!/bin/bash

if [ "$1" == "" ]
then
	filename=parms/parms_paper2.txt
else
	filename=$1
fi

if [ "$2" == "" ]
then
	csvfile=tmp.csv
else
	csvfile=$2
fi


echo `date` Compiling
make release
echo `date` Starting faststi
./faststi -f $filename >$csvfile
echo `date` Starting general analysis
Rscript R/generalAnalysis.R $csvfile
echo #################
echo # PARAMETER FILE: $filename
echo ################
cat $filename
echo `date` All done
