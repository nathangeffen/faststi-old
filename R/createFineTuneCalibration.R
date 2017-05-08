

print(sprintf("Reading file bestParmsInitial.csv"))
inp = read.csv("bestParmsInitial.csv", TRUE)

parameters = ""
for (i in seq(nrow(inp))) {
    parameters = paste(parameters, toString(inp[i,1]), toString(inp[i,2]), "\n")
}
write(parameters, file="tmp1.txt")

######
# Fit better SD_RELATIONSHIP_PERIOD and MEAN_RELATIONSHIP_PERIOD

sd = inp[inp$parameter=="SD_RELATIONSHIP_PERIOD",][1,2]
mn = inp[inp$parameter=="MEAN_RELATIONSHIP_PERIOD",][1,2]
sdseq = seq(sd-0.4, sd+0.4, 0.1)
if (min(sdseq) < 0.1) {
  sdseq = seq(0.1,max(sdseq),0.1)
}
mnseq = seq(mn-0.5,mn+0.4,0.1)

lines = paste("MEAN_RELATIONSHIP_PERIOD", spaceSeparatedList(mnseq), "!")
lines = paste(lines, "\nSD_RELATIONSHIP_PERIOD", spaceSeparatedList(sdseq), "+ !\n")

write(lines, file="tmp2.txt")
system("cat parms/parms_100k_base_calibrate.txt tmp1.txt tmp2.txt >tmp_parms_mn_sd.txt")
system("./faststi -f tmp_parms_mn_sd.txt >tmp5.csv")
system("Rscript R/calibrate.R tmp5.csv >tmp6.csv")
cmd = 'grep CALIBRATE tmp6.csv |  awk \'{print "CALIBRATE,"$3",2017.000,PARAMETER,"}\' | uniq | tail -1 | grep -f - tmp5.csv | awk -F "," \'{print $5" "$6}\' >tmp7.txt'
system(cmd)

######
# Fit better SD_SINGLE_PERIOD and MEAN_SINGLE_PERIOD
sd = inp[inp$parameter=="SD_SINGLE_PERIOD",][1,2]
mn = inp[inp$parameter=="MEAN_SINGLE_PERIOD",][1,2]
sdseq = seq(sd-0.4, sd+0.4, 0.1)
if (min(sdseq) < 0.1) {
  sdseq = seq(0.1,max(sdseq),0.1)
}
mnseq = seq(mn-0.5,mn+0.4,0.1)

lines = paste("MEAN_SINGLE_PERIOD", spaceSeparatedList(mnseq), "!")
lines = paste(lines, "\nSD_SINGLE_PERIOD", spaceSeparatedList(sdseq), "+ !\n")

write(lines, file="tmp2.txt")
system("cat parms/parms_100k_base_calibrate.txt tmp1.txt tmp7.txt tmp2.txt >tmp_parms_mn_sd.txt")
system("./faststi -f tmp_parms_mn_sd.txt >tmp5.csv")
system("Rscript R/calibrate.R tmp5.csv >tmp6.csv")
cmd = 'grep CALIBRATE tmp6.csv |  awk \'{print "CALIBRATE,"$3",2017.000,PARAMETER,"}\' | uniq | tail -1 | grep -f - tmp5.csv | awk -F "," \'{print $5" "$6}\' >tmp8.txt'
system(cmd)

#### 
# Fit better SCALE_SINGLE_PERIOD_INITIAL SCALE_SINGLE_PERIOD_DURING

si = inp[inp$parameter=="SCALE_SINGLE_PERIOD_INITIAL",][1,2]
sd = inp[inp$parameter=="SCALE_SINGLE_PERIOD_DURING",][1,2]

siseq = seq(max(0.01,si-0.1), min(si+0.1, 1.0), 0.01)
sdseq = seq(max(0.01,sd-0.1), min(sd+0.1, 1.0), 0.01)


lines = paste("SCALE_SINGLE_PERIOD_INITIAL", spaceSeparatedList(siseq), "!")
lines = paste(lines, "\nSCALE_SINGLE_PERIOD_DURING", spaceSeparatedList(sdseq), "+ !\n")

write(lines, file="tmp2.txt")
system("cat parms/parms_100k_base_calibrate.txt tmp1.txt tmp7.txt tmp8.txt tmp2.txt >tmp_parms_mn_sd.txt")
system("./faststi -f tmp_parms_mn_sd.txt >tmp5.csv")
system("Rscript R/calibrate.R tmp5.csv >tmp6.csv")
cmd = 'grep CALIBRATE tmp6.csv |  awk \'{print "CALIBRATE,"$3",2017.000,PARAMETER,"}\' | uniq | tail -1 | grep -f - tmp5.csv | awk -F "," \'{print $5" "$6}\' >tmp9.txt'
system(cmd)

#### 
# Fit better SCALE_SINGLE_PERIOD_ZERO_DAYS_INITIAL SCALE_SINGLE_PERIOD_ZERO_DAYS_DURING

si = inp[inp$parameter=="SCALE_SINGLE_PERIOD_ZERO_DAYS_INITIAL",][1,2]
sd = inp[inp$parameter=="SCALE_SINGLE_PERIOD_ZERO_DAYS_DURING",][1,2]

siseq = seq(max(0.01,si-0.1), min(si+0.1, 1.0), 0.01)
sdseq = seq(max(0.01,sd-0.1), min(sd+0.1, 1.0), 0.01)


lines = paste("SCALE_SINGLE_PERIOD_ZERO_DAYS_INITIAL", spaceSeparatedList(siseq), "!")
lines = paste(lines, "\nSCALE_SINGLE_PERIOD_ZERO_DAYS_DURING", spaceSeparatedList(sdseq), "+ !\n")

write(lines, file="tmp2.txt")
system("cat parms/parms_100k_base_calibrate.txt tmp1.txt tmp7.txt tmp8.txt tmp9.txt tmp2.txt >tmp_parms_mn_sd.txt")
system("./faststi -f tmp_parms_mn_sd.txt >tmp5.csv")
system("Rscript R/calibrate.R tmp5.csv >tmp6.csv")
cmd = 'grep CALIBRATE tmp6.csv |  awk \'{print "CALIBRATE,"$3",2017.000,PARAMETER,"}\' | uniq | tail -1 | grep -f - tmp5.csv | awk -F "," \'{print $5" "$6}\' >tmp10.txt'
system(cmd)

#### 
# Fit better SCALE_RELATIONSHIP_PERIOD_INITIAL SCALE_RELATIONSHIP_PERIOD_DURING

si = inp[inp$parameter=="SCALE_RELATIONSHIP_PERIOD_INITIAL",][1,2]
sd = inp[inp$parameter=="SCALE_RELATIONSHIP_PERIOD_DURING",][1,2]

siseq = seq(max(0.01,si-0.1), min(si+0.1, 1.0), 0.01)
sdseq = seq(max(0.01,sd-0.1), min(sd+0.1, 1.0), 0.01)


lines = paste("SCALE_RELATIONSHIP_PERIOD_INITIAL", spaceSeparatedList(siseq), "!")
lines = paste(lines, "\nSCALE_RELATIONSHIP_PERIOD_DURING", spaceSeparatedList(sdseq), "+ !\n")

write(lines, file="tmp2.txt")
system("cat parms/parms_100k_base_calibrate.txt tmp1.txt tmp7.txt tmp8.txt tmp9.txt tmp10.txt tmp2.txt >tmp_parms_mn_sd.txt")
system("./faststi -f tmp_parms_mn_sd.txt >tmp5.csv")
system("Rscript R/calibrate.R tmp5.csv >tmp6.csv")
cmd = 'grep CALIBRATE tmp6.csv |  awk \'{print "CALIBRATE,"$3",2017.000,PARAMETER,"}\' | uniq | tail -1 | grep -f - tmp5.csv | awk -F "," \'{print $5" "$6}\' >tmp11.txt'
system(cmd)

system("cat tmp7.txt tmp8.txt tmp9.txt tmp10.txt tmp11.txt >finalFittedParms.txt")
