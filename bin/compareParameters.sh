./faststi -f parms/parms_test.txt >tmp.csv
head -1 tmp.csv >tmp1.csv
grep ",MEAN_RELATIONSHIP_PERIOD," tmp.csv >> tmp1.csv
grep ",SD_RELATIONSHIP_PERIOD," tmp.csv  >>tmp1.csv
grep ",MEAN_SINGLE_PERIOD," tmp.csv  >>tmp1.csv
grep ",SD_SINGLE_PERIOD," tmp.csv  >>tmp1.csv
grep ",SCALE_RELATIONSHIP_PERIOD_INITIAL," tmp.csv  >>tmp1.csv
grep ",SCALE_RELATIONSHIP_PERIOD_DURING," tmp.csv  >>tmp1.csv
grep ",SCALE_SINGLE_PERIOD_ZERO_DAYS_INITIAL," tmp.csv  >>tmp1.csv
grep ",SCALE_SINGLE_PERIOD_ZERO_DAYS_DURING," tmp.csv  >>tmp1.csv
grep ",MEAN_CASUAL_SEX," tmp.csv  >>tmp1.csv
grep ",SD_CASUAL_SEX," tmp.csv  >>tmp1.csv
grep ",PARTNERSHIPS," tmp.csv | grep 2020.000 >>tmp1.csv
grep ",CASUAL," tmp.csv | grep 2020.000 >>tmp1.csv
grep ",POOR," tmp.csv | grep 2020.000 >>tmp1.csv
grep ",PREVALENCE," tmp.csv | grep 2017.000 | head -1 >>tmp1.csv
grep ",PREVALENCE," tmp.csv | grep 2020.000 >>tmp1.csv
grep ",SINGLES," tmp.csv >>tmp1.csv
grep -E ",MATINGPOOL,|,BREAKUPS," tmp.csv >>tmp1.csv
Rscript R/analyzeParameters.R $1 >> tmp1.csv
cp tmp1.csv $1.csv
