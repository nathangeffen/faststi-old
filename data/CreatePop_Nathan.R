
path <- "C:/Users/Stefan Scholz/Dropbox/2_Arbeit/Planet_Romeo/ContactMatrix/"

pop.size <- 10000 # Set population size 

df <- read.csv(paste(path,"input.csv",sep=""), header=T)
df$AGE.GROUPS <- rep(seq(1,5), each = 4)
df$ACTIVITY <- rep(seq(1,4),5)

df$Nlow <- round(df$lSUPPLY * pop.size, 0)
df$Nrep <- round(df$pSUPPLY * pop.size, 0)

pop.low <- df[rep(row.names(df), df$Nlow), c(2:3)]
pop.low$ID <- seq(1, dim(pop.low)[1],1)
pop.low <- pop.low[,c(3,1,2)]
colnames(pop.low) <- c("ID", "AGE", "ACT")
#pop.low$SEXOR <- 0

write.csv(pop.low, paste(path,"pop_low.csv",sep=""), row.names = FALSE)

pop.rep <- df[rep(row.names(df), df$Nrep), c(2:3)]
pop.rep$ID <- seq(1, dim(pop.rep)[1],1)
pop.rep <- pop.rep[,c(3,1,2)]
colnames(pop.rep) <- c("ID", "AGE", "ACT")
#pop.rep$SEXOR <- 0

write.csv(pop.rep, paste(path,"pop_rep.csv",sep=""), row.names = FALSE)

