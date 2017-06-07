echo `date` Compiling
make release
echo `date` Starting faststi
./faststi -f parms/parms_paper2.txt >tmp.csv
echo `date` Starting general analysis
Rscript R/generalAnalysis.R tmp.csv >results.csv
echo `date` All done
