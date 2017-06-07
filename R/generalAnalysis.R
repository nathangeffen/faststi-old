args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  filename = "tmp.csv"
} else {
  filename = args[1]
}

inp = read.csv(filename, header=TRUE)

prevalences = inp[inp$Desc2=="PREVALENCE",]
analysisDate = max(prevalences$Date)
prevalences = prevalences[prevalences$Date==analysisDate,]

meanPrevalences = aggregate(prevalences$Value, by=list(prevalences$Name),
                    FUN=mean)

partnerships = inp[inp$Desc2=="PARTNERSHIPS",]
partnerships = partnerships[partnerships$Date==analysisDate,]

meanPartnerships = aggregate(partnerships$Value, by=list(partnerships$Name),
                              FUN=mean)

timings = inp[inp$Desc1=="TIMING",]
timings = timings[timings$Date==analysisDate,]

meanTimings = aggregate(timings$Value, by=list(timings$Name),
                             FUN=mean)


sprintf("%s,%f", meanPrevalences$Group.1, meanPrevalences$x)
sprintf("%s,%f", meanPartnerships$Group.1, meanPartnerships$x)
sprintf("%s,%f", meanTimings$Group.1, meanTimings$x)

