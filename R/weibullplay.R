shapes=seq(0.1,10,0.1)
scales=seq(1,100,1)

df = data.frame(shape=0,scale=0,mean=0,median=0,min=0,max=0)

for (i in shapes) {
  for (j in scales) {
    vals = rweibull(100,i, j)
    df <- rbind(df, c(i, j, mean(vals), median(vals), min(vals), max(vals)))
  } 
}


dist <- function(shape, scale) {
  vals = rweibull(100, shape, scale)
  print(mean(vals))
  print(median(vals))
  print(min(vals))
  print(max(vals))
  print(quantile(vals, probs=c(0.25,0.75)))
  hist(vals)
}
