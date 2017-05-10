library(xtable)
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  filename = "tmp.csv"
} else {
  filename = args[1]
}

# Timings
inp = read.csv(filename, TRUE)
timings = inp[inp$Desc1=="TIMING" & inp$Desc2=="AFTER",]
mean_timings = aggregate(timings$Value, by=list(timings$Name),FUN=mean)
colnames(mean_timings) <- c("Simulation","Avg. Time")
xtable(mean_timings)

# Prevalences

confinterval <- function(vals) {
  tryCatch(
    {
      as.double(quantile(vals,probs=c(0.025,0.975)))
      # x = t.test(vals)
      # x$conf.int
    },
    error = function(e) {
      NA
    }
  )
}


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

calcPrevalenceForDate <- function(date_of_analysis) {
  prevalences = inp[inp$Desc1 == "ANALYSIS",]
  prevalences = prevalences[grep("PREVALENCE|MALE",prevalences$Desc2),]
  prevalences = prevalences[prevalences$Date==date_of_analysis,]
  mean_prevalences =  aggregate(prevalences$Value,
                                by=list(prevalences$Name,
                                        prevalences$Desc2),
                                FUN=mean)
  sd_prevalences = aggregate(prevalences$Value,
                             by=list(prevalences$Name,
                                     prevalences$Desc2),
                              FUN=sd)
  ci_prevalences_02_5 = aggregate(prevalences$Value,
                                  by=list(prevalences$Name,
                                          prevalences$Desc2),
                           FUN=confinterval_02_5)
  ci_prevalences_97_5 = aggregate(prevalences$Value,
                                  by=list(prevalences$Name,
                                          prevalences$Desc2),
                               FUN=confinterval_97_5)
  output = data.frame(simulation=mean_prevalences$Group.1,
                    measure=mean_prevalences$Group.2,
                    mean=mean_prevalences$x,
                    sd=sd_prevalences$x,
                    ci_02_5=ci_prevalences_02_5$x,
                    ci_97_5=ci_prevalences_97_5$x)
  output[is.na(output$mean)==FALSE,]

  #output[is.nan(output$mean)==FALSE,]
}

initialprevs = calcPrevalenceForDate(2017)
finalprevs = calcPrevalenceForDate(2020)

format_prevs <- function(prevs) {
  categories = c("ALL", "Males", "Females", "MSM", "WSW",
                 "Male 15-19", "Female 15-19", "Male 20-24", "Female 20-24",
                 "Male 25-29", "Female 25-29", "Male 30-34", "Female 30-34",
                 "Male 35-39", "Female 35-39", "Male 40-44", "Female 40-44",
                 "Male 45-49", "Female 45-49")
  vals = rep(-1, length(categories))
  tex_col_names = c("RPM 0.001", "RKPM 0.001", "CSPM 0.001", "BFPM 0.001", "BLOSSOM 0.001",
                  "RPM 0.01", "RKPM 0.01", "CSPM 0.01", "BFPM 0.01", "BLOSSOM 0.01",
                  "RPM 0.1", "RKPM 0.1", "CSPM 0.1", "BFPM 0.1", "BLOSSOM 0.1")
  r_col_names = gsub(" ","_", tex_col_names)
  r_col_names = gsub("[.]", "_", r_col_names)
  output_table = data.frame(row.names = categories)

  tex_row_names = row.names(output_table)
  s = toupper(tex_row_names)
  s = gsub(" ","_", s)
  #s = gsub("-","-0", s)
  s = gsub("MALE", "MALE_PREVALENCE_AGE", s)
  s = gsub("WSW", "WSW_PREVALENCE", s)
  s = gsub("MSM", "MSM_PREVALENCE", s)
  s = gsub("ALL", "PREVALENCE", s)
  r_row_names = gsub("MALE_AGES", "MALE_PREVALENCE", s)
  for (i in c(0:length(r_row_names))) {
    for (j in c(0:length(r_col_names))) {
      expected = format(round(subset(prevs,measure==r_row_names[i] &
                            simulation==r_col_names[j])$mean[1], 3))
      stddev = format(round(subset(prevs,measure==r_row_names[i] &
                              simulation==r_col_names[j])$sd[1], 3))
      c_025 = format(round(subset(prevs,measure==r_row_names[i] &
                                   simulation==r_col_names[j])$ci_02_5[1], 3))
      c_975 = format(round(subset(prevs,measure==r_row_names[i] &
                                  simulation==r_col_names[j])$ci_97_5[1], 3))
      output_table[tex_row_names[i],  tex_col_names[j] ] = paste(expected, " ",
                                                               "[", c_025, ";", c_975,"]",
                                                               sep="")
    }
  }
  output_table
}

output_table_initial = format_prevs(initialprevs)
xtable(output_table_initial)

output_table_final = format_prevs(finalprevs)
print(xtable(output_table_final))

# Scores

scores = inp[inp$Desc2 == "SCORE",]
scores = scores[scores$Date==2020,]
mean_scores =  aggregate(scores$Value, by=list(scores$Name),
                              FUN=mean)
sd_scores = aggregate(scores$Value, by=list(scores$Name),
                           FUN=sd)

ci_scores_02_5 = aggregate(scores$Value, by=list(scores$Name,
                                                 scores$Desc2),
                                FUN=confinterval_02_5)
ci_scores_97_5 = aggregate(scores$Value, by=list(scores$Name,
                                              scores$Desc2),
                                FUN=confinterval_97_5)
print("SCORES")
print(mean_scores)
print(sd_scores)
print(ci_scores_02_5)
print(ci_scores_97_5)

# Poor

poor = inp[inp$Desc2 == "POOR",]
poor = poor[poor$Date==2020,]
mean_poor =  aggregate(poor$Value, by=list(poor$Name,poor$Desc2),
                         FUN=mean)
sd_poor = aggregate(poor$Value, by=list(poor$Name,poor$Desc2),
                      FUN=sd)
ci_poor_02_5 = aggregate(poor$Value, by=list(poor$Name,poor$Desc2),
                           FUN=confinterval_02_5)
ci_poor_97_5 = aggregate(poor$Value, by=list(poor$Name,poor$Desc2),
                           FUN=confinterval_97_5)

print("POOR")
print(mean_poor)
print(sd_poor)
print(ci_poor_02_5)
print(ci_poor_97_5)


# Histograms to view normality

# cspm_0_01 = prevalences[prevalences$V1=="CSPM_0_01",]
# cspm_0_01 = cspm_0_01[cspm_0_01$V3=="PREVALENCE",]
# hist(cspm_0_01$V6)
#
# rpm_0_01 = prevalences[prevalences$V1=="RPM_0_01",]
# rpm_0_01 = rpm_0_01[rpm_0_01$V3=="PREVALENCE",]
# hist(rpm_0_01$V6)
#
# rpm_0_1 = prevalences[prevalences$V1=="RPM_0_1",]
# rpm_0_1 = rpm_0_1[rpm_0_1$V3=="PREVALENCE",]
# hist(rpm_0_1$V6)
