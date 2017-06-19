#!/bin/bash

if [ "$1" == "" ]
then
	filename=parms/parms_paper2.txt
else
	filename=$1
fi
echo `date` Compiling
make release
echo `date` Starting faststi
./faststi -f $filename >tmp.csv
echo `date` Starting general analysis
Rscript R/generalAnalysis.R tmp.csv
echo `date` All done
