echo `date` Compiling
make release
echo `date` Starting faststi
./faststi -f parms/parms_20k_calibrate.R >tmp.csv
echo `date` Starting grep
grep -E "Name|PARAMETER|MATING|BREAKUP" tmp.csv >tmp2.csv
echo `date` Starting calibrate
Rscript R/calibrate.R tmp2.csv >tmp3.csv
echo `date` Getting best set of parameters
grep CALIBRATE tmp3.csv |  awk '{print "CALIBRATE,"$3",2017.000,PARAMETER,"}' | uniq | tail -1 | grep -f - tmp2.csv | awk -F "," 'BEGIN{print "parameter,value"}{print $5","$6}' >bestParmsInitial.csv
echo `date` Fine tuning parammeters
Rscript R/createFineTuneCalibration.R
echo `date` Prefix parms for 20k runs
cat finalFittedParms.txt parms/parms_20k.txt >tmp_inp_20k.txt
echo `date` Executing 20k runs
./faststi -f tmp_inp_20k.txt >tmp_out_20k.csv
echo `date` Preparing 20k analysis
head -1 tmp_out_20k.csv >tmp20k.csv
grep -E "\,ANALYSIS\,|\,TIMING\," tmp_out_20k.csv >>tmp20k.csv
echo `date` Analysing 20k runs
Rscript R/analyseOutput.R tmp20k.csv >analysis20k.txt
echo `date` Prefix parms for 40m runs
cat finalFittedParms.txt parms/parms_40m.txt >tmp_inp_40m.txt
echo `date` Executing 40m runs
./faststi -f tmp_inp_40m.txt >tmp_out_40m.csv
echo `date` Preparing 40j analysis
head -1 tmp_out_40m.csv >tmp40m.csv
grep -E "\,ANALYSIS\,|\,TIMING\," tmp_out_40m.csv >>tmp40m.csv
echo `date` Analysing 40m runs
Rscript R/analyseOutput.R tmp40m.csv >analysis40m.txt
echo `date` All done
