#loading required libraries 

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()
source("http://bioconductor.org/biocLite.R")
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38"))

BiocManager::install(c("RMariaDB"))

BiocManager::install(c("GenomicFeatures"))
BiocManager::install(c("VariantAnnotation"))

library(GenomicFeatures)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationDbi)
library(parallel)



# Load Human Genome Data
genome <- BSgenome.Hsapiens.UCSC.hg38
#download.file("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz","gencode.v39.primary_assembly.annotation.gtf.gz")

# Load Gencode v26 Annotations
txdb <- makeTxDbFromGFF('data/gencode.v26.primary_assembly.annotation.gtf.gz')
# Save as SQLite Database (if needed)
saveDb(txdb, 'data/txdb.gencode26.sqlite')


# Load Gencode Database
txdb <- loadDb('data/txdb.gencode26.sqlite')
ensgene <- txdb
chr.list <- paste('chr',1:22,sep='')
ensgene.sub <- keepSeqlevels(ensgene,chr.list)
ensgene <- ensgene.sub

cds_seqs <- extractTranscriptSeqs(Hsapiens, cdsBy(ensgene, by="tx",use.names=TRUE))
utr.grange <-threeUTRsByTranscript(ensgene, use.names=TRUE)
fiveutr.grange<-fiveUTRsByTranscript(ensgene, use.names=TRUE)
introns.grange<- intronsByTranscript(ensgene, use.names=TRUE)
cds.grange<- cdsBy(ensgene, by="tx",use.names=TRUE)


###get DNAstring of UTR seqs
five_utr_seqs <- extractTranscriptSeqs(Hsapiens, utr.grange)
three_utr_seqs <- extractTranscriptSeqs(Hsapiens, utr.grange)
five_utr_seqs <- extractTranscriptSeqs(Hsapiens, fiveutr.grange)





keys <- names(cds_seqs)
cols <- columns(ensgene)
all.df <- select(ensgene, keys = keys, columns = cols, keytype="TXNAME")
## this has length of 304803 so most contain NAs or duplicates?
##cds_seqs length is 28856
all.df.sub <- all.df[which(all.df$TXNAME!='NA'),]
###remove dulpicate entries?
#all.df.sub <- all.df[unique(all.df$TXNAME),]
#u.keys<-unique(names(cds_seqs))
#u.all.df <- select(ccds, keys = keys, columns=cols,keytype="TXNAME")



seqs_uORF <- function(seqs,x) {
  print(x)
  a <- unlist(seqs[x])
  name <- names(seqs[x])
  ind <- unlist(gregexpr('ATG',as.character(a)))
  #augs.ind <- ind[which(ind%%3==1)]
  ind.1 <- unlist(gregexpr('TAG|TAA|TGA',as.character(a)))
  #stop.ind <- ind.1[which(ind.1%%3==1)]
  

  if (length(ind)>0 & length(ind.1)>0 & any((ind.1-ind-3)%%3==0) & any(ind.1-ind>0)) {
    paste('there is a utr_uORF')
  } else 
    paste('not any utr_uORF')

}

### for getting the content for every 10 bp

get_GC_content <- function(seqs,x){
  print(x)
  a <- unlist(seqs[x])
  
  if (length(a)<10){
    return(NA)
  }
  starts <- seq(1, length(a), by = 10)
  n <- length(starts)    # Find the length of the vector "starts"
  result<-c()
  for (i in 1:n) {
   chunk <- subseq(a, i, i+9)
    content <- alphabetFrequency(chunk, baseOnly=TRUE, as.prob = TRUE)
    GC<- content[3]+content[2]
    result[i]<-(GC)
}
 return(result)
  
}


get_AU_content <- function(seqs,x){
  print(x)
  a <- unlist(seqs[x])
  
  if (length(a)<10){
    return(NA)
  }
  starts <- seq(1, length(a), by = 10)
  n <- length(starts)    # Find the length of the vector "starts"
  result<-c()
  for (i in 1:n) {
   chunk <- subseq(a, i, i+9)
    content <- alphabetFrequency(chunk, baseOnly=TRUE, as.prob = TRUE)
    AU<- content[1]+content[4]
    result[i]<-(AU)
}
 return(result)
  
}

get_UC_content <- function(seqs,x){
  print(x)
  a <- unlist(seqs[x])
  
  if (length(a)<10){
    return(NA)
  }
  starts <- seq(1, length(a), by = 10)
  n <- length(starts)    # Find the length of the vector "starts"
  result<-c()
  for (i in 1:n) {
   chunk <- subseq(a, i, i+9)
    content <- alphabetFrequency(chunk, baseOnly=TRUE, as.prob = TRUE)
    UC<- content[4]+content[2]
    result[i]<-(UC)
}
 return(result)
  
}


### for the getting the content for every sequence

get_GC_content <- function(seqs,x){
  print(x)
  a <- unlist(seqs[x])
  content <- alphabetFrequency(a, baseOnly=TRUE, as.prob = TRUE)
  GC<- content[3]+content[2]
  
  return(GC)
  
}

get_AU_content <- function(seqs,x){
  print(x)
  a <- unlist(seqs[x])
  content <- alphabetFrequency(a, baseOnly=TRUE, as.prob = TRUE)
  AU<- content[1]+content[4]
  #result[i]<-(AU)
  return(AU)
  
}


get_UC_content <- function(seqs,x){
  print(x)
  a <- unlist(seqs[x])
  content <- alphabetFrequency(a, baseOnly=TRUE, as.prob = TRUE)
  UC<- content[4]+content[2]
  #result[i]<-(UC)
  return(UC)
  
}

fiveutrseqs.uORF <- mclapply(1:length(five_utr_seqs), function(x)
{
     print(x)
     seqs_uORF(five_utr_seqs,x)

},mc.cores=8)


threeutrseqs.uORF <- mclapply(1:length(three_utr_seqs), function(x)
{
     print(x)
     seqs_uORF(three_utr_seqs,x)

},mc.cores=8)

threeUTR_GC_content<- mclapply(1:length(three_utr_seqs), function(x)
{
     print(x)
     get_GC_content(three_utr_seqs,x)

},mc.cores=8)


threeUTR_AU_content<- mclapply(1:length(three_utr_seqs), function(x)
{
     print(x)
     get_AU_content(three_utr_seqs,x)

},mc.cores=4)

threeUTR_UC_content<- mclapply(1:length(three_utr_seqs), function(x)
{
     print(x)
     get_UC_content(three_utr_seqs,x)

},mc.cores=4)


### for cdsseqs GC content

cdsseqs_GC_content_10bp<- mclapply(1:length(cds_seqs), function(x)
{
     print(x)
     get_GC_content(cds_seqs,x)

},mc.cores=1)



fiveUTR_GC_content_10bp<- mclapply(1:length(five_utr_seqs), function(x)
{
     print(x)
     get_GC_content(five_utr_seqs,x)

},mc.cores=2)

fiveUTR_AU_content_10bp<- mclapply(1:length(five_utr_seqs), function(x)
{
     print(x)
     get_AU_content(five_utr_seqs,x)

},mc.cores=2)




#amount of windows with AC>=.8

get_windows<- function(seqs,x){
  
  a<-seqs[[x]]
  paste(length(a[a>=1]))
}



threeUTR_100_AU<- mclapply(1:length(three_utr_seqs), function(x)
{
     get_windows(threeUTR_AU_content,x)

},mc.cores=1)


threeUTR_100_GC<- mclapply(1:length(three_utr_seqs), function(x)
{
     get_windows(threeUTR_GC_content,x)

},mc.cores=1)

threeUTR_100_UC<- mclapply(1:length(three_utr_seqs), function(x)
{
     get_windows(threeUTR_UC_content,x)

},mc.cores=1)

cdsseqs_GC_content <- cdsseqs_GC_content_10bp

cdsseqs_100_GC<- mclapply(1:length(cds_seqs), function(x)
{
     get_windows(cdsseqs_GC_content,x)

},mc.cores=1)





get_AU_last200content <- function(x){
  
  print (x)
  a <- unlist(three_utr_seqs[x])
  
  if (length(a)<200){
    return(NA)}
  
  start<-length(a)-200+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[1]+content[4]
  return(AU)
  
}
get_GC_last200content <- function(x){
    
    print (x)
    a <- unlist(three_utr_seqs[x])
    
    if (length(a)<200){
        return(NA)}
    
    
    start<-length(a)-200+1
    end<-length(a)
    b<-subseq(a, start=start, end=end)
    content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
    AU<-content[3]+content[2]
    return(AU)
    
}

get_UC_last200content <- function(x){
    
    print (x)
    a <- unlist(three_utr_seqs[x])
    
    if (length(a)<200){
        return(NA)}
   
    
    start<-length(a)-200+1
    end<-length(a)
    b<-subseq(a, start=start, end=end)
    content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
    AU<-content[4]+content[2]
    return(AU)
    
}


get_AU_first200content <- function(x){
    
    print (x)
    a <- unlist(three_utr_seqs[x])
    
    if (length(a)<200){
        return(NA)}
    
    #start<-length(a)-200
    #end<-length(a)
    start <- 1
    end <- 200
    b<-subseq(a, start=start, end=end)
    content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
    AU<-content[1]+content[4]
    return(AU)
    
}

get_GC_first200content <- function(x){
    
    print (x)
    a <- unlist(three_utr_seqs[x])
    
    if (length(a)<200){
        return(NA)}

    #start<-length(a)-200
    #end<-length(a)
    start <- 1
    end <- 200
    b<-subseq(a, start=start, end=end)
    content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
    AU<-content[3]+content[2]
    return(AU)
    
}

get_UC_first200content <- function(x){
    
    print (x)
    a <- unlist(three_utr_seqs[x])
    
    if (length(a)<200){
        return(NA)}
  
    #start<-length(a)-200
    #end<-length(a)
    start <- 1
    end <- 200
    b<-subseq(a, start=start, end=end)
    content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
    AU<-content[4]+content[2]
    return(AU)
    
}




ThreeUTR_AUcontentlast200<- unlist(lapply(1:length(three_utr_seqs), get_AU_last200content))
ThreeUTR_GCcontentlast200<- unlist(lapply(1:length(three_utr_seqs), get_GC_last200content))
ThreeUTR_UCcontentlast200<- unlist(lapply(1:length(three_utr_seqs), get_UC_last200content))

ThreeUTR_AUcontentfirst200<- unlist(lapply(1:length(three_utr_seqs), get_AU_first200content))
ThreeUTR_GCcontentfirst200<- unlist(lapply(1:length(three_utr_seqs), get_GC_first200content))
ThreeUTR_UCcontentfirst200<- unlist(lapply(1:length(three_utr_seqs), get_UC_first200content))






#get_cds_length <- function(x){
  #print(x)
  #a <- unlist(cds_seqs[x])
 #l <- length(a)
  
  #return(l)
#}
#cds_length<- mclapply(1:length(cds_seqs), get_cds_length,mc.cores=8)

get_utr_length <- function(x){
  
  print(x)
  a <- unlist(three_utr_seqs[x])
  l <- length(a)
  
  return(l)
}
threeutr_length<- mclapply(1:length(three_utr_seqs), get_utr_length,mc.cores=8)





  

