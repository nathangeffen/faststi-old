library("xtable")

confinterval_02_5 <- function(vals) {
  tryCatch(
    {
      confinterval(vals)[1]
    },
    error = function(e) {
      NA
    }
  )
}

confinterval_97_5 <- function(vals) {
  tryCatch(
    {
      confinterval(vals)[2]
    },
    error = function(e) {
      NA
    }
  )
}



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


ci_prevalences_02_5 = aggregate(prevalences$Value,
                                by=list(prevalences$Name,
                                        prevalences$Desc2),
                                FUN=confinterval_02_5)
ci_prevalences_97_5 = aggregate(prevalences$Value,
                                by=list(prevalences$Name,
                                        prevalences$Desc2),
                                FUN=confinterval_97_5)


print("Mean prevalence")
sprintf("%s,%f", meanPrevalences$Group.1, meanPrevalences$x)
print("0.025 CI prevalence")
sprintf("%s,%f", ci_prevalences_02_5$Group.1, ci_prevalences_02_5$x)
print("0.0975 CI prevalence")
sprintf("%s,%f", ci_prevalences_97_5$Group.1, ci_prevalences_97_5$x)
print("Mean partnerships")
sprintf("%s,%f", meanPartnerships$Group.1, meanPartnerships$x)
print("Mean timings")
sprintf("%s,%f", meanTimings$Group.1, meanTimings$x)
xtable(meanPrevalences)
xtable(ci_prevalences_02_5)
xtable(ci_prevalences_97_5)
xtable(meanPartnerships)
xtable(meanTimings)
