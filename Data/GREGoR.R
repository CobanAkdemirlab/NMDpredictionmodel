library(GenomicRanges)
gr_all <- do.call(c, rds.merged)
gr_all
length(gr_all)
names(gr_all) <- make.unique(names(gr_all))
gregor <- as.data.frame(gr_all)
#6127
