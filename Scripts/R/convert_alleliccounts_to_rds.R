args <- commandArgs(trailing = TRUE)

IN  = args[1]
OUT = args[2]

x = read.delim(IN, comment.char = "@", stringsAsFactors = F)
saveRDS(x, file = OUT)
