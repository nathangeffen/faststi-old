# Finds sets of parameters that produce stable matings and breakups.

getsingles <- function(index)
{
    singles[singles$Num==index,]
}


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  filename = "tmp.csv"
} else {
  filename = args[1]
}

print(sprintf("Reading file: %s", filename))
inp = read.csv(filename, TRUE)
print("Getting parameters")
parameters = inp[inp$Desc1=="PARAMETER",]
print("Getting singles")
singles = inp[inp$Desc2 == "SINGLES" & inp$Date >= 2017.0,]

# sd_matingpool =  aggregate(matingpool$Value, by=list(matingpool$Num), FUN=sd)
# sim_num = which(sd_matingpool$x==min(sd_matingpool$x)) - 1
# best_parameters_mating_sd = parameters[parameters$Num==sim_num,]
# 
# sd_breakups = aggregate(breakups$Value, by=list(breakups$Num), FUN=sd)
# sim_num = which(sd_breakups$x==min(sd_breakups$x)) - 1
# best_parameters_breakups_sd = parameters[parameters$Num==sim_num,]
# 
# sd_combined = sd_matingpool + sd_breakups
# sim_num = which(sd_combined$x==min(sd_combined$x)) - 1
# best_parameters_combined_sd = parameters[parameters$Num==sim_num,]

print("Beginning tests of linearality")
lm_singles = by(singles,singles$Num, 
                   function(x) {
                     # print(sprintf("Singles num: %d",x$Num[1]))
                     abs(lm(x$Value~x$Date,data=x)$coef[2])
                    } )
sim_num = as.double(toString(which(lm_singles==min(lm_singles)) 
                             - 1))
best_parameters_singles_lm = parameters[parameters$Num==sim_num,]
print("Best singles (lm)")
print(lm_singles[sim_num+1])
print(best_parameters_singles_lm)

print("Beginning tests of max-min")
mm_singles = by(singles,singles$Num, 
                function(x) {
                  max(x$Value) - min(x$value)
                } )
sim_num = as.double(toString(which(mm_singles==min(mm_singles)) 
                             - 1))
best_parameters_singles_mm = parameters[parameters$Num==sim_num,]
print("Best singles (max - min)")
print(mm_singles[sim_num+1])
print(best_parameters_singles_mm)
