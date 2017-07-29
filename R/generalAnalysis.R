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

analyzeVar <- function(varName, analysisDate, description, descField, input)
{
  if(missing(analysisDate)) {analysisDate = 1}
  if(missing(description)) {description = varName}
  if(missing(descField)) {descField = 2}
  if(missing(input)) {input = inp}
  if (descField == 1) {
    values = input[input$Desc1==varName,]
  } else {
    values = input[input$Desc2==varName,]
  }
  if (nrow(values) == 0) {return(0)}
  if (analysisDate == 1) {
    analysisDate = max(values$Date)
  } else {
    analysisDate = min(values$Date)
  }
  values = values[values$Date==analysisDate,]
  meanValues = aggregate(values$Value,
                         by=list(values$Name),
                         FUN=mean)
  ci_02_5 = aggregate(values$Value,
                      by=list(values$Name),
                      FUN=confinterval_02_5)
  ci_97_5 = aggregate(values$Value,
                      by=list(values$Name),
                      FUN=confinterval_97_5)

  #newrow = c(meanValues$Group.1[1], description[1],
  #           as.double(meanValues$x[1]), as.double(ci_02_5$x[1]), as.double(ci_97_5$x[1]))
  #print(newrow)
  #df <- rbind(df, newrow)

  print(sprintf("%s,%s,%.3f,%.3f,%.3f", meanValues$Group.1, description, meanValues$x,
          ci_02_5$x,ci_97_5$x), quote=FALSE)
}



args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  filename = "tmp.csv"
} else {
  filename = args[1]
}

inp = read.csv(filename, header=TRUE)
inp = inp[inp$Desc1 != "PARAMETER",]
inp$Value = as.double(as.character(inp$Value))

analyzeVar("PREVALENCE", 0, "INITIAL PREVALENCE")
analyzeVar("PREVALENCE", description="FINAL PREVALENCE")
analyzeVar("AGE_RANGE_PREVALENCE", 0, "INITIAL PREVALENCE TRUNCATED AGE")
analyzeVar("AGE_RANGE_PREVALENCE", 1, "FINAL PREVALENCE TRUNCATED AGE")
analyzeVar("FEMALE_PREVALENCE_AGE_25-29", 0, "INITIAL FEMALE_PREVALENCE_AGE_25-29")
analyzeVar("FEMALE_PREVALENCE_AGE_15-29", description="FINAL FEMALE_PREVALENCE_AGE_15-29")
analyzeVar("FEMALE_PREVALENCE_AGE_20-24", description="FINAL FEMALE_PREVALENCE_AGE_20-24")
analyzeVar("FEMALE_PREVALENCE_AGE_25-29", description="FINAL FEMALE_PREVALENCE_AGE_25-29")
analyzeVar("FEMALE_PREVALENCE_AGE_30-34", description="FINAL FEMALE_PREVALENCE_AGE_30-34")
analyzeVar("FEMALE_PREVALENCE_AGE_35-39", description="FINAL FEMALE_PREVALENCE_AGE_35-39")
analyzeVar("FEMALE_PREVALENCE_AGE_40-44", description="FINAL FEMALE_PREVALENCE_AGE_40-44")
analyzeVar("FEMALE_PREVALENCE_AGE_45-49", description="FINAL FEMALE_PREVALENCE_AGE_45-49")
analyzeVar("MALE_PREVALENCE_AGE_15-29", description="FINAL MALE_PREVALENCE_AGE_15-29")
analyzeVar("MALE_PREVALENCE_AGE_20-24", description="FINAL MALE_PREVALENCE_AGE_20-24")
analyzeVar("MALE_PREVALENCE_AGE_25-29", description="FINAL MALE_PREVALENCE_AGE_25-29")
analyzeVar("MALE_PREVALENCE_AGE_30-34", description="FINAL MALE_PREVALENCE_AGE_30-34")
analyzeVar("MALE_PREVALENCE_AGE_35-39", description="FINAL MALE_PREVALENCE_AGE_35-39")
analyzeVar("MALE_PREVALENCE_AGE_40-44", description="FINAL MALE_PREVALENCE_AGE_40-44")
analyzeVar("MALE_PREVALENCE_AGE_45-49", description="FINAL MALE_PREVALENCE_AGE_45-49")

analyzeVar("MSM_PREVALENCE",0,"INITIAL MSM PREVALENCE")
analyzeVar("MSM_PREVALENCE", description="FINAL MSM PREVALENCE")
analyzeVar("WSW_PREVALENCE",0,"INITIAL WSW PREVALENCE")
analyzeVar("WSW_PREVALENCE", description="FINAL WSW PREVALENCE")
analyzeVar("PARTNERSHIPS", 0, "INITIAL PARTNERSHIPS")
analyzeVar("PARTNERSHIPS", description="FINAL PARTNERSHIPS")
analyzeVar("CASUAL")
analyzeVar("POOR")
analyzeVar("SCORE")
analyzeVar("SINGLES", 0, "INITIAL SINGLES")
analyzeVar("SINGLES", description="FINAL SINGLES")
analyzeVar("TIMING", descField = 1)
