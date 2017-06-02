args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  splits = 10
} else {
  splits = as.integer(args[1])
}

if (length(args) < 2) {
  filename = "tmp.csv"
} else {
  filename = args[2]
}

print(sprintf("Reading file: %s", filename))
inp = read.csv(filename, TRUE)

last = max(inp$Num)

splitSize = as.integer(last / splits)
if (splitSize * splits < last) {
  splits = splits + 1
}

for (i in seq(splits)) {
    from = (i - 1) * splitSize
    to = i * splitSize
    out = inp[inp$Num >= from & inp$Num < to, ]
    outfile = paste("out_", toString(i), ".csv", sep="")
    write.csv(out, file=outfile)
}
