library(ape)

tr <- read.tree("out.best.dnd")

utr <- unroot(tr)

write.tree(utr, "out.best.unrooted.dnd")