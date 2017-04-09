# Example R script to interact with Faststi (C++ microsimulation of STI epidemics) 
# Faststi takes a parameter file as input and outputs a CSV file with analysis of the epidemic 

# General format of program
# condition = false
# while condition is not met
#   write file for input for Faststi
#   execute Faststi
#   analyse the csv file that Faststi produced

# Input and output files for Faststi 
filenameToFaststi = "inputToFaststi.csv"
filenameFromFaststi = "outputFromFaststi.csv"

# Function to write Faststi's input parameter file
writeFaststiInput <- function(name, forceInfectionFrom, forceInfectionTo, step)
{
  lines = ""
  # In case you want faststi to append to its CSV file on each execution, then
  # you need a unique name for each simulation so that you can aggregate by simulation run 
  # However, if faststi overwrites its output CSV file everytime, then the next line can be 
  # simplified to: 
  # simulation_name = name
  simulation_name = paste(name, toString(forceInfectionFrom), sep="")
  
  
  # Check the method setDefaultParameters in the file src/parameter.hh to 
  # see what all the parameters are that Faststi can take.
  line = sprintf("SIMULATION_NAME %s\n", simulation_name)
  
  lines = paste(line, "NUM_AGENTS 20000\n")
  
  # You might only want to run the simulation once, or some other number of times, in which 
  # case adjust this parameter
  lines = paste(line, "NUM_SIMULATIONS 100\n")
  
  # I have four cores on my machine, so 4 threads makes sense. Adjust as you see fit.
  lines = paste(lines, "NUM_THREADS 4\n")
  lines = paste(lines, "START_DATE 2017\nEND_DATE 2020\nMATCH_NEIGHBORS 30\n")
  lines = paste(lines, "MATCH_CLUSTERS 200\nMATCH_EVENT CSPM\n")
  
  # We produce output every 110 days. Adjust to your needs.
  lines = paste(lines, "ANALYZE_DURING_SIM 110\n")
  
  # In this toy example, we fiddle with the infectiousness rates.
  # But in a more sophisticated example, e.g. to calibrate the model, we'd 
  # have our adjustable parameters here.
  
  # Now this is quite sophisticated and you might not need it.
  # If you specify three values for some parameters and you put a "@" at the line
  # each iteration of the simulation will vary the value across the range of the three values
  # where the first value is the beginning of the range, the second value is the end of the
  # range (and not included in the range) and the third value is the amount to step
  # by for each simulation. Think of it as being like Python's range function.
  # So in the examples below, each of these parameters will increment on each iteration
  # of the simulation. If it reaches the highest value of the range before the last simulation
  # it wraps to the beginning of the range.
  line = sprintf("HET_MALE_INFECTIOUSNESS %.5f %.5f %.5f @\n", 
                 forceInfectionFrom, forceInfectionTo, step)
  lines = paste(lines, line)
  line = sprintf("HET_FEMALE_INFECTIOUSNESS %.5f %.5f %.5f @\n", 
                 forceInfectionFrom / 2.0, forceInfectionTo / 2.0, step / 2.0)
  lines = paste(lines, line)
  line = sprintf("HOM_MALE_INFECTIOUSNESS %.5f %.5f %.5f @\n", 
                 forceInfectionFrom * 2.0, forceInfectionTo * 2.0, step * 2.0)
  lines = paste(lines, line)  
  line = sprintf("HOM_FEMALE_INFECTIOUSNESS %.5f %.5f %.5f @\n", 
                 forceInfectionFrom / 2.0, forceInfectionTo / 2.0, step / 2.0)
  lines = paste(lines, line)  
  write(lines, file=filenameToFaststi)
}

# This is a toy example just to demonstrate. For calibration more sophisticated 
# initialization is probably needed.
prevalence = 0
infectiousnessFrom = 0.001
step = 0.0001
infectiousnessTo = infectiousnessFrom + 100 * step

# This while loop shows the general format. We execute until Faststi 
# produces a set of parameters we're happy with. To avoid an endless loop, it might be 
# a good idea to 

while(prevalence < 0.3) {
  # Write the parameter file
  writeFaststiInput("RTEST_", infectiousnessFrom, infectiousnessTo, step)
  # Prepare and the execute the Linux shell command to execute Faststi 
  cmd = sprintf("./bin/faststi -f %s >%s", filenameToFaststi, filenameFromFaststi)
  system(cmd)
  
  # Read the output of Faststi. This step can be slow because 100 executions of 
  # the simulation produces a lot of output.
  inp = read.csv(filenameFromFaststi, FALSE)
  
  # Filter out everything we don't need
  outputPrevalences = inp[inp$V3 == "PREVALENCE" & inp$V5 == 2020,]
  
  # Find the maximum prevalence
  prevalence = max(outputPrevalences$V6)
  # Find the row in the data frame with the highest prevalence
  row_num = which(outputPrevalences$V6==prevalence)
  
  # (This awful next line of code is the only way I could find to get a scalar
  # from a data frame. Feel free to suggest a more stylish way.)
  # Each row has the simulation number. We need the simulation number to find out
  # what the parameter was that it was executed with. bestSimulation gets that simulation
  # number.
  bestSimulation = as.integer(toString(outputPrevalences[row_num,4])) 
  
  # Now we find the value of the parameter we passed.
  bestParameter =  inp[inp$V4 == "HET_MALE_INFECTIOUSNESS" & inp$V2 == bestSimulation,]$V5
  # Print our info
  print(sprintf("Prevalence: %.4f Infectiousness: %.4f", prevalence, bestParameter))
  
  # Update our parameters
  infectiousnessFrom = infectiousnessTo
  infectiousnessTo = infectiousnessFrom + 100 * step
}
