#
setwd('/path/to/TOR.csv')

files <- list.files('.','.vcf')
het_get <- data.frame()
for (i in 1: length(files)){
print(i)
TOR <- read.table(files[i],header=TRUE)
ind <- which(TOR[,ncol(TOR)]=='0/1')
TOR.sub <- TOR[ind,]
colnames(TOR.sub)[ncol(TOR.sub)]<-'sample'
het_get <- rbind(het_get, TOR.sub)
}
