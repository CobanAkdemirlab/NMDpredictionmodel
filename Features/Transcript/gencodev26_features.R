#libraries
library(GenomicFeatures)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationDbi)
library(parallel)
library(Biostrings)
library(RMariaDB)
library(GenomicRanges)
genome <- BSgenome.Hsapiens.UCSC.hg38

####
txdb = makeTxDbFromGFF('~/data/gencode.v26.primary_assembly.annotation.gtf.gz')
# optional loading if dataset is available
# saveDb(txdb, '~/txdb.gencode26.sqlite')
#start of the sequencing features extraction
ensgene <- txdb
chr.list <- paste('chr',1:22,sep='')
ensgene.sub <- keepSeqlevels(ensgene,chr.list)
ensgene <- ensgene.sub

cds_seqs <- extractTranscriptSeqs(Hsapiens, cdsBy(ensgene, by="tx",use.names=TRUE))
utr.grange <-threeUTRsByTranscript(ensgene, use.names=TRUE)
fiveutr.grange<-fiveUTRsByTranscript(ensgene, use.names=TRUE)
introns.grange<- intronsByTranscript(ensgene, use.names=TRUE)
cds.grange<- cdsBy(ensgene, by="tx",use.names=TRUE)
exons.grange<-exonsBy(ensgene, use.names=TRUE)
###get DNAstring of UTR seqs
cds_seqs <- extractTranscriptSeqs(Hsapiens,
                                  cdsBy(txdb, by="tx", use.names=TRUE))
prot_seqs <- translate(cds_seqs)

# Raw exon-level sequences (DNAStringSetList)
threeutr_parts <- getSeq(BSgenome.Hsapiens.UCSC.hg38, utr.grange)
fiveutr_parts <- getSeq(BSgenome.Hsapiens.UCSC.hg38, fiveutr.grange)

# Collapse exon sequences per transcript with reduce() + xscat
three_utr_seqs <- unlist(endoapply(threeutr_parts, function(x) do.call(xscat, as.list(x))))
five_utr_seqs  <- unlist(endoapply(fiveutr_parts, function(x) do.call(xscat, as.list(x))))

# Keep transcript names
names(three_utr_seqs) <- names(utr.grange)
names(five_utr_seqs)  <- names(fiveutr.grange)

keys <- names(cds_seqs)
cols <- columns(ensgene)
all.df <- select(ensgene, keys = keys, columns = cols, keytype="TXNAME")
all.df.sub <- all.df[which(all.df$TXNAME!='NA'),]

## get the length of 3'UTR
get_utr_length <- function(x){
  
  print(x)
  a <- unlist(three_utr_seqs[x])
  l <- length(a)
  
  return(l)
}
library(parallel)
threeutr_length<- mclapply(1:length(three_utr_seqs), get_utr_length,mc.cores=8)

### get the length of 5'utr
get_5utr_length <- function(x){
  
  print(x)
  a <- unlist(five_utr_seqs[x])
  l <- length(a)
  
  return(l)
}
fiveutr_length<- mclapply(1:length(five_utr_seqs), get_5utr_length,mc.cores=8)

### is there any uORF in a sequences
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
### for the getting the content for every sequence
## GC content
get_GC_content <- function(seqs,x){
  print(x)
  a <- unlist(seqs[x])
  content <- alphabetFrequency(a, baseOnly=TRUE, as.prob = TRUE)
  GC<- content[3]+content[2]
  
  return(GC)
  
}

## AU content

get_AU_content <- function(seqs,x){
  print(x)
  a <- unlist(seqs[x])
  content <- alphabetFrequency(a, baseOnly=TRUE, as.prob = TRUE)
  AU<- content[1]+content[4]
  #result[i]<-(AU)
  return(AU)
  
}
## UC content

get_UC_content <- function(seqs,x){
  print(x)
  a <- unlist(seqs[x])
  content <- alphabetFrequency(a, baseOnly=TRUE, as.prob = TRUE)
  UC<- content[4]+content[2]
  #result[i]<-(UC)
  return(UC)
  
}






### UORF in fiveutrseqs
library(parallel)

fiveutrseqs.uORF <- mclapply(1:length(five_utr_seqs), function(x)
{
  print(x)
  seqs_uORF(five_utr_seqs,x)
  
},mc.cores=8)


### UORF in threeutrseqs # no need--

threeutrseqs.uORF <- mclapply(1:length(three_utr_seqs), function(x)
{
  print(x)
  seqs_uORF(three_utr_seqs,x)
  
},mc.cores=8)

## GC content in three UTR
threeUTR_GC_content<- mclapply(1:length(three_utr_seqs), function(x)
{
  print(x)
  get_GC_content(three_utr_seqs,x)
  
},mc.cores=8)

## AU content in three UTR

threeUTR_AU_content<- mclapply(1:length(three_utr_seqs), function(x)
{
  print(x)
  get_AU_content(three_utr_seqs,x)
  
},mc.cores=4)

## UC content in three UTR
threeUTR_UC_content<- mclapply(1:length(three_utr_seqs), function(x)
{
  print(x)
  get_UC_content(three_utr_seqs,x)
  
},mc.cores=4)

## GC content in five UTR
fiveUTR_GC_content<- mclapply(1:length(five_utr_seqs), function(x)
{
  print(x)
  get_GC_content(five_utr_seqs,x)
  
},mc.cores=8)

## AU content in five UTR

fiveUTR_AU_content<- mclapply(1:length(five_utr_seqs), function(x)
{
  print(x)
  get_AU_content(five_utr_seqs,x)
  
},mc.cores=4)

## ## UC content in five UTR

fiveUTR_UC_content<- mclapply(1:length(five_utr_seqs), function(x)
{
  print(x)
  get_UC_content(five_utr_seqs,x)
  
},mc.cores=4)

### for cdsseqs GC content

cdsseqs_GC_content<- mclapply(1:length(cds_seqs), function(x)
{
  print(x)
  get_GC_content(cds_seqs,x)
  
},mc.cores=1)

cdsseqs_AU_content<- mclapply(1:length(cds_seqs), function(x)
{
  print(x)
  get_AU_content(cds_seqs,x)
  
},mc.cores=1)

cdsseqs_UC_content<- mclapply(1:length(cds_seqs), function(x)
{
  print(x)
  get_UC_content(cds_seqs,x)
  
},mc.cores=1)


######################## Three UTR first/last 100bp ###################
get_AU_last100content <- function(x){
  
  print (x)
  a <- unlist(three_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  start<-length(a)-100+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[1]+content[4]
  return(AU)
  
}
get_GC_last100content <- function(x){
  
  print (x)
  a <- unlist(three_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  
  start<-length(a)-100+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[3]+content[2]
  return(AU)
  
}
get_UC_last100content <- function(x){
  
  print (x)
  a <- unlist(three_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  
  start<-length(a)-100+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[4]+content[2]
  return(AU)
  
}

get_AU_first100content <- function(x){
  
  print (x)
  a <- unlist(three_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  #start<-length(a)-200
  #end<-length(a)
  start <- 1
  end <- 100
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[1]+content[4]
  return(AU)
  
}
get_GC_first100content <- function(x){
  
  print (x)
  a <- unlist(three_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  #start<-length(a)-200
  #end<-length(a)
  start <- 1
  end <- 100
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[3]+content[2]
  return(AU)
  
}

get_UC_first100content <- function(x){
  
  print (x)
  a <- unlist(three_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  #start<-length(a)-200
  #end<-length(a)
  start <- 1
  end <- 100
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[4]+content[2]
  return(AU)
  
}

ThreeUTR_AUcontentlast100<- unlist(lapply(1:length(three_utr_seqs), get_AU_last100content))
ThreeUTR_GCcontentlast100<- unlist(lapply(1:length(three_utr_seqs), get_GC_last100content))
ThreeUTR_UCcontentlast100<- unlist(lapply(1:length(three_utr_seqs), get_UC_last100content))

ThreeUTR_AUcontentfirst100<- unlist(lapply(1:length(three_utr_seqs), get_AU_first100content))
ThreeUTR_GCcontentfirst100<- unlist(lapply(1:length(three_utr_seqs), get_GC_first100content))
ThreeUTR_UCcontentfirst100<- unlist(lapply(1:length(three_utr_seqs), get_UC_first100content))


#######################200

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

## fiveutr content
########################threeutr100#####
get_AU_5utrlast100content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  start<-length(a)-100+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[1]+content[4]
  return(AU)
  
}
get_GC_5utrlast100content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  
  start<-length(a)-100+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[3]+content[2]
  return(AU)
  
}
get_UC_5utrlast100content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  
  start<-length(a)-100+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[4]+content[2]
  return(AU)
  
}

get_AU_5utrfirst100content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  #start<-length(a)-200
  #end<-length(a)
  start <- 1
  end <- 100
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[1]+content[4]
  return(AU)
  
}
get_GC_5utrfirst100content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  #start<-length(a)-200
  #end<-length(a)
  start <- 1
  end <- 100
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[3]+content[2]
  return(AU)
  
}

get_UC_5utrfirst100content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  #start<-length(a)-200
  #end<-length(a)
  start <- 1
  end <- 100
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[4]+content[2]
  return(AU)
  
}

fiveUTR_AUcontentlast100<- unlist(lapply(1:length(five_utr_seqs), get_AU_5utrlast100content))
fiveUTR_GCcontentlast100<- unlist(lapply(1:length(five_utr_seqs), get_GC_5utrlast100content))
fiveUTR_UCcontentlast100<- unlist(lapply(1:length(five_utr_seqs), get_UC_5utrlast100content))

fiveUTR_AUcontentfirst100<- unlist(lapply(1:length(five_utr_seqs), get_AU_5utrfirst100content))
fiveUTR_GCcontentfirst100<- unlist(lapply(1:length(five_utr_seqs), get_GC_5utrfirst100content))
fiveUTR_UCcontentfirst100<- unlist(lapply(1:length(five_utr_seqs), get_UC_5utrfirst100content))


# maybe 200

get_AU_5utrlast200content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
  if (length(a)<200){
    return(NA)}
  
  start<-length(a)-200+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[1]+content[4]
  return(AU)
  
}
get_GC_5utrlast200content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
  if (length(a)<200){
    return(NA)}
  
  
  start<-length(a)-200+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[3]+content[2]
  return(AU)
  
}

get_UC_5utrlast200content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
  if (length(a)<200){
    return(NA)}
  
  
  start<-length(a)-200+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[4]+content[2]
  return(AU)
  
}


get_AU_5utrfirst200content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
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

get_GC_5utrfirst200content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
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

get_UC_5utrfirst200content <- function(x){
  
  print (x)
  a <- unlist(five_utr_seqs[x])
  
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




FiveUTR_AUcontentlast200<- unlist(lapply(1:length(five_utr_seqs), get_AU_5utrlast200content))
FiveUTR_GCcontentlast200<- unlist(lapply(1:length(five_utr_seqs), get_GC_5utrlast200content))
FiveUTR_UCcontentlast200<- unlist(lapply(1:length(five_utr_seqs), get_UC_5utrlast200content))

FiveUTR_AUcontentfirst200<- unlist(lapply(1:length(five_utr_seqs), get_AU_5utrfirst200content))
FiveUTR_GCcontentfirst200<- unlist(lapply(1:length(five_utr_seqs), get_GC_5utrfirst200content))
FiveUTR_UCcontentfirst200<- unlist(lapply(1:length(five_utr_seqs), get_UC_5utrfirst200content))

## cds



get_cds_length <- function(x){
print(x)
  a <- unlist(cds_seqs[x])
  l <- length(a)

return(l)
}
cds_length<- mclapply(1:length(cds_seqs), get_cds_length,mc.cores=8)



cdsseqs_GC_content<- mclapply(1:length(cds_seqs), function(x)
{
  print(x)
  get_GC_content(cds_seqs,x)
  
},mc.cores=1)

cdsseqs_AU_content<- mclapply(1:length(cds_seqs), function(x)
{
  print(x)
  get_AU_content(cds_seqs,x)
  
},mc.cores=1)

cdsseqs_UC_content<- mclapply(1:length(cds_seqs), function(x)
{
  print(x)
  get_UC_content(cds_seqs,x)
  
},mc.cores=1)


#
## GC content 

get_AU_cdslast100content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  start<-length(a)-100+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[1]+content[4]
  return(AU)
  
}
get_GC_cdslast100content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  
  start<-length(a)-100+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[3]+content[2]
  return(AU)
  
}
get_UC_cdslast100content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  
  start<-length(a)-100+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[4]+content[2]
  return(AU)
  
}

get_AU_cdsfirst100content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  #start<-length(a)-200
  #end<-length(a)
  start <- 1
  end <- 100
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[1]+content[4]
  return(AU)
  
}
get_GC_cdsfirst100content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  #start<-length(a)-200
  #end<-length(a)
  start <- 1
  end <- 100
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[3]+content[2]
  return(AU)
  
}

get_UC_cdsfirst100content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
  if (length(a)<100){
    return(NA)}
  
  #start<-length(a)-200
  #end<-length(a)
  start <- 1
  end <- 100
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[4]+content[2]
  return(AU)
  
}

cdsseq_AUcontentlast100<- unlist(lapply(1:length(cds_seqs), get_AU_cdslast100content))
cdsseq_GCcontentlast100<- unlist(lapply(1:length(cds_seqs), get_GC_cdslast100content))
cdsseq_UCcontentlast100<- unlist(lapply(1:length(cds_seqs), get_UC_cdslast100content))

cdsseq_AUcontentfirst100<- unlist(lapply(1:length(cds_seqs), get_AU_cdsfirst100content))
cdsseq_GCcontentfirst100<- unlist(lapply(1:length(cds_seqs), get_GC_cdsfirst100content))
cdsseq_UCcontentfirst100<- unlist(lapply(1:length(cds_seqs), get_UC_cdsfirst100content))


#200
get_AU_cdslast200content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
  if (length(a)<200){
    return(NA)}
  
  start<-length(a)-200+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[1]+content[4]
  return(AU)
  
}
get_GC_cdslast200content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
  if (length(a)<200){
    return(NA)}
  
  
  start<-length(a)-200+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[3]+content[2]
  return(AU)
  
}

get_UC_cdslast200content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
  if (length(a)<200){
    return(NA)}
  
  
  start<-length(a)-200+1
  end<-length(a)
  b<-subseq(a, start=start, end=end)
  content <- alphabetFrequency(b, baseOnly=TRUE, as.prob = TRUE)
  AU<-content[4]+content[2]
  return(AU)
  
}


get_AU_cdsfirst200content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
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

get_GC_cdsfirst200content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
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

get_UC_cdsfirst200content <- function(x){
  
  print (x)
  a <- unlist(cds_seqs[x])
  
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




cdsseq_AUcontentlast200<- unlist(lapply(1:length(cds_seqs), get_AU_cdslast200content))
cdsseq_GCcontentlast200<- unlist(lapply(1:length(cds_seqs), get_GC_cdslast200content))
cdsseq_UCcontentlast200<- unlist(lapply(1:length(cds_seqs), get_UC_cdslast200content))

cdsseq_AUcontentfirst200<- unlist(lapply(1:length(cds_seqs), get_AU_cdsfirst200content))
cdsseq_GCcontentfirst200<- unlist(lapply(1:length(cds_seqs), get_GC_cdsfirst200content))
cdsseq_UCcontentfirst200<- unlist(lapply(1:length(cds_seqs), get_UC_cdsfirst200content))

# 
####number of exons
cds_exons <- function(x)
{
  print(x)
  name <- names(cds_seqs[x])
  df <- all.df.sub[which(all.df.sub$TXNAME==name),]
  cds.df<-df[!is.na(df$CDSID),]
  num.exons<-nrow(cds.df)
  cds_exon_size<-(cds.df$CDSEND)-(cds.df$CDSSTART)+1
  names(cds_exon_size)<-c(1:num.exons)
  exons<-cumsum(cds_exon_size)
  paste(exons, collapse=',')
  
}

cds.exons<-mclapply(1:length(cds_seqs), cds_exons,mc.cores=8)


cds_length <- function(x)
{
  print(x)
  name <- names(cds_seqs[x])
  df <- all.df.sub[which(all.df.sub$TXNAME==name),]
  cds.df<-df[!is.na(df$CDSID),]
  cds_exon_size<-(cds.df$CDSEND)-(cds.df$CDSSTART)+1
  s<-sum(cds_exon_size)
  return(s)
}

cds.length<- mclapply(1:length(cds_seqs), cds_length,mc.cores=1)


NC_length <- function(x)
{
  print(x)
  name <- names(cds_seqs[x])
  df <- all.df.sub[which(all.df.sub$TXNAME==name),]
  cds.df<-df[is.na(df$CDSID),]
  cds_exon_size<-(cds.df$EXONEND)-(cds.df$EXONSTART)+1
  s<-sum(cds_exon_size)
  return(s)
}
noncoding.length<- mclapply(1:length(cds_seqs), NC_length,mc.cores=1)

#hall.df <- select(refgene, keys = keys, columns=cols, keytype="TXNAME")
## this has length of 304803 so most contain NAs or duplicates?
##cds_seqs length is 28856
all.df.sub <- all.df[which(all.df$TXNAME!='NA'),]

exonnum <- function(x)
{
  print(x)
  name <- names(cds_seqs[x])
  df <- all.df[which(all.df$TXNAME==name),]
  cds.df<-df[!is.na(df$CDSID),]
  num.exons<-nrow(cds.df)
  return(num.exons)
}
exon_count<-mclapply(1:length(cds_seqs), exonnum,mc.cores=8)

NC_exons<- function(x){
  print (x)
  name <- names(cds_seqs[x])
  df <- all.df[which(all.df$TXNAME==name),]
  
  NCcds.df<-df[is.na(df$CDSID),]
  NCnum.exons<-nrow(NCcds.df)
  return(NCnum.exons)
  
}
NCexonsnum<-mclapply(1:length(cds_seqs), NC_exons,mc.cores=8)




### what is the location of initiation alternative codons

downstream_start <- function(x)
{
  print(x)
  a <- unlist(cds_seqs[x])
  name <- names(cds_seqs[x])
  ind <- unlist(gregexpr('ATG',as.character(a)))
  plus1.ind <- ind[which(ind%%3==1)]
  if (length(plus1.ind)>0) {
    startcodons <- plus1.ind
    paste(plus1.ind,collapse = ",")
  } else 
    paste(NA)
}

downstream.start <- mclapply(1:length(cds_seqs), downstream_start,mc.cores=8)


#
hits <- unique(queryHits(findOverlaps(utr.grange,introns.grange)))
txnames.3utr.introns <- names(utr.grange) [hits]

# if there is any intron in 5utr
hits <- unique(queryHits(findOverlaps(fiveutr.grange,introns.grange)))
txnames.5utr.introns <- names(fiveutr.grange) [hits]

# if there is any intron in 
hits <- unique(queryHits(findOverlaps(cds.grange,introns.grange)))
txnames.cdsseq.introns <- names(cds.grange) [hits]

### organizing the files for extraction

cds.seq.fr <- data.frame(txnames=names(cds_seqs), cds_length=unlist(cds_length), cds_exons=unlist(cds.exons), 
                         noncoding.length=unlist(noncoding.length), exon_count=unlist(exon_count),NCexonsnum=unlist(NCexonsnum), 
                         downstream_start=unlist(downstream.start), cdsseqs_GC_content=unlist(cdsseqs_GC_content),
                         cdsseqs_UC_content=unlist(cdsseqs_UC_content), cdsseqs_AU_content=unlist(cdsseqs_AU_content), 
                         cdsseq_AUcontentlast100=cdsseq_AUcontentlast100, cdsseq_UCcontentlast100=cdsseq_UCcontentlast100, 
                         cdsseq_GCcontentlast100=cdsseq_GCcontentlast100, cdsseq_AUcontentfirst100=cdsseq_AUcontentfirst100, 
                         cdsseq_UCcontentfirst100=cdsseq_UCcontentfirst100, cdsseq_GCcontentfirst100=cdsseq_GCcontentfirst100,
                         cdsseq_AUcontentlast200=cdsseq_AUcontentlast200, cdsseq_UCcontentlast200=cdsseq_UCcontentlast200, 
                         cdsseq_GCcontentlast200=cdsseq_GCcontentlast200, cdsseq_AUcontentfirst200=cdsseq_AUcontentfirst200, 
                         cdsseq_UCcontentfirst200=cdsseq_UCcontentfirst200, cdsseq_GCcontentfirst200=cdsseq_GCcontentfirst200)

three.utr.fr <- data.frame(txnames=names(three_utr_seqs), threeutrseqs.uORF=unlist(threeutrseqs.uORF), 
                           threeUTR_AU_content=unlist(threeUTR_AU_content), threeUTR_GC_content=unlist(threeUTR_GC_content),
                           threeUTR_UC_content=unlist(threeUTR_UC_content), threeUTR_length=unlist(threeutr_length),
                           ThreeUTR_AUcontentlast200=ThreeUTR_AUcontentlast200,ThreeUTR_AUcontentfirst200=ThreeUTR_AUcontentfirst200,
                           ThreeUTR_GCcontentfirst200=ThreeUTR_GCcontentfirst200, ThreeUTR_GCcontentlast200=ThreeUTR_GCcontentlast200,
                           ThreeUTR_UCcontentfirst200=ThreeUTR_UCcontentfirst200,ThreeUTR_UCcontentlast200=ThreeUTR_UCcontentlast200,
                           ThreeUTR_UCcontentlast100=ThreeUTR_UCcontentlast100, ThreeUTR_UCcontentfirst100=ThreeUTR_UCcontentfirst100,
                           ThreeUTR_GCcontentfirst100=ThreeUTR_GCcontentfirst100, ThreeUTR_GCcontentlast100=ThreeUTR_GCcontentlast100, 
                           ThreeUTR_AUcontentfirst100=ThreeUTR_AUcontentfirst100, ThreeUTR_AUcontentlast100=ThreeUTR_AUcontentlast100)


five.utr.fr <- data.frame(txnames=names(five_utr_seqs),fiveutrseqs.uORF=unlist(fiveutrseqs.uORF), 
                          fiveutr_length=unlist(fiveutr_length), fiveUTR_AU_content=unlist(fiveUTR_AU_content), 
                          fiveUTR_GC_content=unlist(fiveUTR_GC_content), fiveUTR_UC_content=unlist(fiveUTR_UC_content),
                          FiveUTR_AUcontentlast200=FiveUTR_AUcontentlast200, FiveUTR_AUcontentfirst200=FiveUTR_AUcontentfirst200,
                          FiveUTR_UCcontentlast200=FiveUTR_UCcontentlast200, FiveUTR_UCcontentfirst200=FiveUTR_UCcontentfirst200,
                          FiveUTR_GCcontentlast200=FiveUTR_GCcontentlast200, FiveUTR_GCcontentfirst200=FiveUTR_GCcontentfirst200,
                          fiveUTR_AUcontentlast100=fiveUTR_AUcontentlast100, fiveUTR_AUcontentfirst100= fiveUTR_AUcontentfirst100,
                          fiveUTR_UCcontentlast100=fiveUTR_UCcontentlast100, fiveUTR_UCcontentfirst100=fiveUTR_UCcontentfirst100,
                          fiveUTR_GCcontentlast100=fiveUTR_GCcontentlast100, fiveUTR_GCcontentfirst100=fiveUTR_GCcontentfirst100)


### first merge cds.seq.fr and three.utr.fr

merged.1 <- merge(cds.seq.fr,three.utr.fr,all.x=TRUE,by='txnames')
merged.2 <- merge(merged.1,five.utr.fr,all.x=TRUE,by='txnames')

merged.2$threeUTR.introns <- rep('NA',nrow(merged.2))
merged.2[which(merged.2$txnames%in%txnames.3utr.introns),'threeUTR.introns'] <- 'There is a 3UTR intron'

merged.2$fiveUTR.introns <- rep('NA',nrow(merged.2))
merged.2[which(merged.2$txnames%in%txnames.5utr.introns),'fiveUTR.introns'] <- 'There is a 5UTR intron'

### get the PTBP1 bindings sites in first 200-400 bp of 3'utr 
utr.resized=resize(utr.grange,width=400,fix='start')
### get the PTBP1 binding sites
PTBP1.1 <- read.table('~/ENCFF907HNN.bed.txt',sep='\t')
PTBP1.2 <- read.table('~/ENCFF130PWU.bed.txt',sep='\t')
PTBP1.3 <- read.table('~/ENCFF100OEX.bed.txt',sep='\t')


PTBP1 <- rbind(PTBP1.1,PTBP1.2,PTBP1.3)


gr2 <- GRanges(seqnames = PTBP1$V1, ranges = IRanges(start = PTBP1$V2, end=PTBP1$V3))

hits <- unique(queryHits(findOverlaps(utr.resized,gr2)))
txnames.3utr.PTBP1 <- names(utr.grange) [hits]
merged.2$threeUTR.PTBP1 <- rep('NA',nrow(merged.2))
merged.2[which(merged.2$txnames%in%txnames.3utr.PTBP1),'threeUTR.PTBP1'] <- 'There is a PTBP1 binding motif'
merged.2$threeUTR.PTBP1 <-
  merged.2$txnames %in% txnames.3utr.PTBP1

save(merged.2, file='hg38_seqfeatures_Gencodev26_2026.RData')

#Sanity check:
dim(merged.2)
#[1] 91907    61
