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
meanPrevalences
meanPartnerships