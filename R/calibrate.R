# Finds sets of parameters that produce stable matings and breakups.

getmatings <- function(index)
{
    matingpool[matingpool$Num==index,]
}

getbreakups <- function(index)
{
  breakups[breakups$Num==index,]
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
print("Getting matings")
matingpool = inp[inp$Desc1 == "MATINGPOOL" & inp$Date >= 2017.0,]
print("Getting breakups")
breakups = inp[inp$Desc1 == "BREAKUPS" & inp$Date >= 2017.0,]

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
lm_matingpool = by(matingpool,matingpool$Num, 
                   function(x) {
                     print(sprintf("Mating num: %d",x$Num[1]))
                     abs(lm(x$Value~x$Date,data=x)$coef[2])
                    } )
sim_num = as.double(toString(which(lm_matingpool==min(lm_matingpool)) 
                             - 1))
best_parameters_mating_lm = parameters[parameters$Num==sim_num,]
print("Best mating")
print(lm_matingpool[sim_num+1])
print(best_parameters_mating_lm)

lm_breakups = by(breakups,breakups$Num, 
                 function(x) {
                   print(sprintf("Breakup num: %d",x$Num[1]))
                   abs(lm(x$Value~x$Date,data=x)$coef[2])
                  } )
sim_num = which(lm_breakups==min(lm_breakups)) - 1
best_parameters_breakups_lm = parameters[parameters$Num==sim_num,]
print("Best breakups")
print(lm_breakups[sim_num+1])
print(best_parameters_breakups_lm)

lm_combined = lm_breakups + lm_matingpool
sim_num = which(lm_combined==min(lm_combined)) - 1
best_parameters_combined_lm = parameters[parameters$Num==sim_num,]
print("Best combined")
print(lm_combined[sim_num+1])
print(best_parameters_combined_lm)

# best_parameters_combined_lm
# best_parameters_mating_lm
# best_parameters_breakups_lm
# best_parameters_combined_sd
# best_parameters_mating_sd
# best_parameters_breakups_sd
# 


