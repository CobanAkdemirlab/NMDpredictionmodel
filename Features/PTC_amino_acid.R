#After running gencodeV26 or gencodeV45 code
#or re-run again
if (!require("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("txdbmaker")

library(ggplot2)
library(gridExtra)
library(plyr)
library(dplyr)
library(reshape2)
library(parallel)
library(plotly)
library (data.table)
library (ggpubr)
library(GenomicFeatures)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationDbi)

genome <- BSgenome.Hsapiens.UCSC.hg38
library(txdbmaker)

txdb <- makeTxDbFromGFF("~/gencode/gencode.v26.primary_assembly.annotation.gtf.gz")


#saveDb(txdb, '~/location/txdb.gencode45.sqlite')
ensgene <- txdb


chr.list <- paste('chr',1:22,sep='')
ensgene.sub <- keepSeqlevels(ensgene,chr.list)
ensgene <- ensgene.sub

cds_seqs <- extractTranscriptSeqs(Hsapiens, cdsBy(ensgene, by="tx",use.names=TRUE))
utr.grange <-threeUTRsByTranscript(ensgene, use.names=TRUE)
fiveutr.grange<-fiveUTRsByTranscript(ensgene, use.names=TRUE)
introns.grange<- intronsByTranscript(ensgene, use.names=TRUE)
exons.grange<-exonsBy(ensgene, use.names=TRUE)

cds.grange<- cdsBy(ensgene, by="tx",use.names=TRUE)


###get DNAstring of UTR seqs
three_utr_seqs <- extractTranscriptSeqs(Hsapiens, utr.grange)
five_utr_seqs <- extractTranscriptSeqs(Hsapiens, fiveutr.grange)
cds_seqs <- extractTranscriptSeqs(Hsapiens,
                                  cdsBy(txdb, by="tx", use.names=TRUE))
prot_seqs <- translate(cds_seqs)


#categorization
#PTC (amino acid)
df$minus1_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.minus1.aa=='G')
df$minus1_is_G[ind.1] <- rep('Glycine',length(ind.1))


df$minus2_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.minus2.aa=='G')
df$minus2_is_G[ind.1] <- rep('Glycine',length(ind.1))

df$plus1_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.plus1.aa=='G')
df$plus1_is_G[ind.1] <- rep('Glycine',length(ind.1))

df$plus2_is_G <- rep('Not Glycine',nrow(df))
ind.1 <- which(df$PTC.plus2.aa=='G')
df$plus2_is_G[ind.1] <- rep('Glycine',length(ind.1))


