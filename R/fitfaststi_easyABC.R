# Example R script to interact with Faststi (C++ microsimulation of STI epidemics) 
# Faststi takes a parameter file as input and outputs a CSV file with analysis of the epidemic 

# General format of program
# condition = false
# while condition is not met
#   write file for input for Faststi
#   execute Faststi
#   analyse the csv file that Faststi produced

# set working directory for relative paths
setwd("~/dev/faststi/")

initial_prevs <- function(msw = 0.01,
                          msm = 0.01,
                          wsm = 0.01,
                          wsw = 0.01) {
  return(as.list(environment()))
}

# Age categories

infection_parms <- function(inf_msw = 0.1,
                            inf_msm = 0.2,
                            inf_wsm = 0.05,
                            inf_wsw = 0.05) {
  return(as.list(environment()))
}


# MIN_AGE_TRUNCATE
# MAX_AGE_TRUNCATE -> can be left out then 12-100
# Sex survey data problematic after 40 or 50
# parameter to ask for analysis only between MIN and MAX
# additional parameter -> follow truncated at start people to the end
# INC_TRUNCATED_AGE -> increment truncated age as simulation progresses
# analyzeTruncatedAge _Init _During _After (default = FALSE)

write_input <- function(name = "RTEST_", 
                        n_simulations = 100, 
                        n_agents = 10000, 
                        start_date = 1980, 
                        end_date = 2020,
                        analysis_step_days = 365,
                        n_threads = 4,
                        prevs = initial_prevs(),
                        parms = infection_parms()) {
  # Write initial rates to file
  with(prevs, {
    initial_rates <- data.frame(AGE = c(14,19,24,29,34,39,44,49,100), # 14 - everyone up to 14.99999, 19 = 15-19.9999
                                MSW = c(rep(msw,8),0),
                                MSM = c(rep(msw,8),0),
                                WSM = c(rep(wsm,8),0),
                                WSW = c(rep(wsw,8),0))
    write.csv(initial_rates, "InitialRates.csv", row.names = FALSE, quote = FALSE)
  })
  # Write parameters to file
  with(parms, {
    lines <- ""
    simulation_name <- paste0(name, toString(inf_msw))
    line <- sprintf("SIMULATION_NAME %s\n", simulation_name)
    lines <- paste0(lines, line)
    lines <- paste0(lines, "NUM_SIMULATIONS ", n_simulations, "\n")
    lines <- paste0(lines, "NUM_AGENTS ", n_agents, "\n")
    lines <- paste0(lines, "NUM_THREADS ", n_threads, "\n")
    lines <- paste0(lines, "START_DATE ", start_date, "\nEND_DATE ", end_date, "\n")
    lines <- paste0(lines, "MATCH_EVENT CSPM\n") # CSPM / RKPM / BLOSSOMV
    lines <- paste0(lines, "MATCH_NEIGHBORS 10\n") # 
    lines <- paste0(lines, "MATCH_CLUSTERS 5\n") #
    lines <- paste0(lines, "ANALYZE_DURING_SIM ", analysis_step_days, "\n")
    lines <- paste0(lines, "HET_MALE_INFECTIOUSNESS ", inf_msw, "\n")
    lines <- paste0(lines, "HOM_MALE_INFECTIOUSNESS ", inf_msm, "\n")
    lines <- paste0(lines, "HET_FEMALE_INFECTIOUSNESS ", inf_wsm, "\n")
    lines <- paste0(lines, "HOM_FEMALE_INFECTIOUSNESS ", inf_wsw, "\n")
    lines <- paste0(lines, "BREAKUP_EVENT LIMIT\nMATING_POOL_EVENT LIMIT\n")
    lines <- paste0(lines, "INITIAL_INFECTION_RATES_CSV InitialRates.csv")
    write(lines, file="FaststiInput.csv")
  })
}

# Test manually with a while loop

mean_final_prevalence <- 0
i <- 1
msw_transprob <- 0.1
step <- 0.01

while(mean_final_prevalence < 0.3) {
  print(paste("Simulation run",i))
  write_input(n_simulations = 10,
              start_date = 1985,
              end_date = 2020,
              prevs = initial_prevs(msw = 0.001, 
                                    msm = 0.005, 
                                    wsm = 0.001, 
                                    wsw = 0),
              parms = infection_parms(inf_msw = msw_transprob,
                                      inf_msm = msw_transprob*2,
                                      inf_wsm = msw_transprob/2,
                                      inf_wsw = msw_transprob/4))
  cmd = sprintf("./bin/faststi -f %s >%s", "FaststiInput.csv", "FaststiOutput.csv")
  system(cmd)
  results <- read.csv("FaststiOutput.csv", TRUE)
  final_prevalences <- dplyr::filter(results, Desc2 == "PREVALENCE", Date == max(Date))
  mean_final_prevalence = mean(final_prevalences$Value)
  print(paste0("MSW infectivity: ",msw_transprob))
  print(paste0("Mean final prevalence: ", mean_final_prevalence))
  
  # Increment counter and set infectivity higher
  i <- i + 1
  msw_transprob <- msw_transprob + step
}

# BREAKUP_EVENT -> FREQUENCY / LIMIT / RANDOM
# MATING_POOL -> FREQUENCY / LIMIT / RANDOM
# above must be the same

set_priors <- function(msw.tprob.mean = 0.1,
                       sd = 0.1,
                       ratio.msmmsw = 2,
                       ratio.wsmmsw = 1/2,
                       ratio.wswmsw = 1/4) {
  priors <- list(prior.msw.tprob = c("lognormal",msw.tprob.mean,sd),
                 prior.msm.tprob = c("lognormal",msw.tprob.mean*ratio.msmmsw,sd),
                 prior.wsm.tprob = c("lognormal",msw.tprob.mean*ratio.wsmmsw,sd),
                 prior.wsw.tprob = c("lognormal",msw.tprob.mean*ratio.wswmsw,sd))
  return(priors)
}

faststi_simulation <- function(seed = 123, param.set) {
  # This is an ugly workaround for the fact that EasyABC only allows parameters
  # for which priors are set
  fixed <- get("fixed_params",env=globalenv())
  
  with(fixed, {
    write_input(n_simulations = n_simulations,
                n_agents = 10000,
                start_date = start_date,
                end_date = end_date,
                prevs = prevs,
                parms = infection_parms(inf_msw = param.set[1],
                                        inf_msm = param.set[2],
                                        inf_wsm = param.set[3],
                                        inf_wsw = param.set[4]))
  })
  time <- toString(lubridate::now())
  outfilename <- paste0("FaststiOutput_",gsub(" ","_",gsub(":","-",time)),".csv")
  cmd = sprintf("./bin/faststi -s %i -f %s >%s", seed, "FaststiInput.csv", outfilename)
  system(cmd)
  results <- read.csv(outfilename, TRUE)
  final_prevalences <- dplyr::filter(results, Desc2 == "PREVALENCE", Date == max(Date))
  mean_final_prevalence = mean(final_prevalences$Value)
  sd_final_prevalence = sd(final_prevalences$Value)
  print(mean_final_prevalence)
  return(c(mean_final_prevalence))
}



fixed_params <-  list(n_simulations = 4,
                      n_threads = 4,
                      start_date = 2015,
                      end_date = 2020,
                      prevs = initial_prevs(msw = 0.01, 
                                            msm = 0.01, 
                                            wsm = 0.01, 
                                            wsw = 0))
priors <- list(prior.msw.tprob = c("unif",0,1),
               prior.msm.tprob = c("unif",0,1),
               prior.wsm.tprob = c("unif",0,1),
               prior.wsw.tprob = c("unif",0,1))

faststi_simulation(param.set = c(0.05,0.05,0.05,0.05))


#priors <- set_priors()
library(EasyABC)
sim.results <- ABC_sequential(model=faststi_simulation,
                              method = "Lenormand",
                             prior=priors,
                             prior_test = "X2 > X1 && X3 < X1 && X3 > X4",
                             summary_stat_target = 0.05,
                             progress_bar=TRUE,
                             nb_simul = 20,
                             alpha = 0.5,
                             p_acc_min = 0.05)
sim.df <- as.data.frame(sim.results)
View(sim.df)

