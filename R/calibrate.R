getmatings <- function(index)
{
    matingpool[matingpool$V4==index,]
}

getbreakups <- function(index)
{
  breakups[breakups$V4==index,]
}


inp = read.csv("tmp.csv", FALSE)
parameters = inp[grep("PARAMETER",inp$V2),]
matingpool = inp[inp$V2=="MATINGPOOL" & inp$V5>=2017.0,]
breakups = inp[inp$V2=="BREAKUPS" & inp$V5>=2017.0,]

sd_matingpool =  aggregate(matingpool$V6, by=list(matingpool$V4), FUN=sd)
sim_num = which(sd_matingpool$x==min(sd_matingpool$x)) - 1
best_parameters_mating_sd = parameters[parameters$V4==sim_num,]

sd_breakups = aggregate(breakups$V6, by=list(breakups$V4), FUN=sd)
sim_num = which(sd_breakups$x==min(sd_breakups$x)) - 1
best_parameters_breakups_sd = parameters[parameters$V4==sim_num,]

sd_combined = sd_matingpool + sd_breakups
sim_num = which(sd_combined$x==min(sd_combined$x)) - 1
best_parameters_combined_sd = parameters[parameters$V4==sim_num,]

lm_matingpool = by(matingpool,matingpool$V4, function(x) {abs(lm(x$V6~x$V5,data=x)$coef[2])} )
sim_num = which(lm_matingpool==min(lm_matingpool)) - 1
best_parameters_mating_lm = parameters[parameters$V4==sim_num,]

lm_breakups = by(breakups,breakups$V4, function(x) {abs(lm(x$V6~x$V5,data=x)$coef[2])} )
sim_num = which(lm_breakups==min(lm_breakups)) - 1
best_parameters_breakups_lm = parameters[parameters$V4==sim_num,]

lm_combined = lm_breakups + lm_matingpool
sim_num = which(lm_combined==min(lm_combined)) - 1
best_parameters_combined_lm = parameters[parameters$V4==sim_num,]

best_parameters_combined_lm[1,4]
best_parameters_mating_lm[1,4]
best_parameters_breakups_lm[1,4]
best_parameters_combined_sd[1,4]
best_parameters_mating_sd[1,4]
best_parameters_breakups_sd[1,4]



