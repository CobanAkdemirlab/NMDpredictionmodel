#library
library(readxl)
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

library(biomaRt)
library(AnnotationDbi)
library(RColorBrewer)


library(ggsignif)
library (mltools)

library(hrbrthemes)
library(viridis)
library(ggh4x)
library(wesanderson)
library(rstatix)
library(varImp)
library(tidyverse)
##

#opitional if re-loading:
#variants.features.fr <- df
#new data

ind <- which(is.na(variants.features.fr$cds_exons))
variants.features.fr <- variants.features.fr[-ind,]
variants.features.fr$coding.pos <-  sapply(1:nrow(variants.features.fr),function(x)
  
{
  print(x)
  ens.annotated <- variants.features.fr$V3[x]
  pos <- strsplit(ens.annotated,':')
  cod.pos <- pos[[1]][grep(variants.features.fr$txnames[x],pos[[1]])+2]
  if(length(grep('c.',cod.pos))==1) {
    if (length(grep('_',cod.pos))==1) {
      cod.pos <- strsplit(cod.pos,'_')[[1]][1]
      as.numeric(gsub("\\D", "", cod.pos))
    } else {
      
      as.numeric(gsub("\\D", "", cod.pos))
    }
  } else 
  {
    'NA'
    
  }
  
  
  
})

variants.features.fr$mut.exon <- sapply(1:nrow(variants.features.fr), function(x) {

  if (!is.na(variants.features.fr$coding.pos[x]) &&
      !is.na(variants.features.fr$cds_exons[x])) {
    
    exon.cor <- as.numeric(unlist(strsplit(
      as.character(variants.features.fr$cds_exons[x]), ',')))
    
    if (variants.features.fr$coding.pos[x] <= exon.cor[1]) {
      return(1)
      
    } else if (length(exon.cor) >= 2) {
      
      for (i in 2:length(exon.cor)) {
        if (variants.features.fr$coding.pos[x] >= exon.cor[i-1] &&
            variants.features.fr$coding.pos[x] <= exon.cor[i]) {
          return(i)
        }
      }
    }
  }

  return(NA_real_)  # ALWAYS numeric
})

variants.features.fr$GENE_ID <-  sapply(1:nrow(variants.features.fr),function(x)
  
{
  
  ens.annotated <- variants.features.fr$V3[x]
  pos <- strsplit(ens.annotated,':')
  gene <- pos[[1]][1]
  
})

##mutated exon length
variants.features.fr$length.mutated.exon <- sapply(1:nrow(variants.features.fr), function(x) {

  if (!is.na(variants.features.fr$cds_exons[x]) &&
      !is.na(variants.features.fr$mut.exon[x])) {

    exon.cor <- as.numeric(unlist(strsplit(
      as.character(variants.features.fr$cds_exons[x]), ',')))

    mut.exon <- variants.features.fr$mut.exon[x]

    if (mut.exon >= 2) {
      return(exon.cor[mut.exon] - exon.cor[mut.exon - 1])
    } else {
      return(exon.cor[mut.exon])
    }
  }

  return(NA_real_)   
})


variants.features.fr$first.100 <- sapply(seq_len(nrow(variants.features.fr)), function(x) {
  val <- as.numeric(variants.features.fr$coding.pos[x])
  if (!is.na(val) && val <= 100) {
    "first 100"
  } else {
    "not first 100"
  }
})

#First 200bp
variants.features.fr$first.200 <- sapply(seq_len(nrow(variants.features.fr)), function(x) {
  val <- as.numeric(variants.features.fr$coding.pos[x])
  if (!is.na(val) && val <= 200) {
    "first 200"
  } else {
    "not first 200"
  }
})


#### long exon 
 variants.features.fr$long.exon <- sapply(1:nrow(variants.features.fr), function(x) {

  if (!is.na(variants.features.fr$cds_exons[x]) &&
      !is.na(variants.features.fr$mut.exon[x])) {

    exon.cor <- as.numeric(unlist(strsplit(
      as.character(variants.features.fr$cds_exons[x]), ',')))

    mut.exon <- variants.features.fr$mut.exon[x]

    if (mut.exon >= 2) {
      length.exon <- exon.cor[mut.exon] - exon.cor[mut.exon - 1]
    } else {
      length.exon <- exon.cor[mut.exon]
    }

    if (length.exon >= 400) {
      return("long.exon")
    } else {
      return("not long.exon")
    }
  }

  return(NA_character_)  # 🔥 critical
})
### PT2 to EJC distance

variants.features.fr$PTC.2.EJC<-  sapply(1:nrow(variants.features.fr),function(x)
  
{
  print(x)
  
  if(!is.na(variants.features.fr$cds_exons[x]) & !variants.features.fr$mut.exon[x]=='NA' & !is.null(variants.features.fr$mut.exon[x][[1]])) {
    exon.cor <- as.numeric(unlist(strsplit(as.character(variants.features.fr$cds_exons[x]),',')))
    mut.exon <- variants.features.fr$mut.exon[x]
    if(mut.exon[[1]] >=2){
      length.exon=exon.cor[mut.exon[[1]]]-as.numeric(variants.features.fr$coding.pos[x])
      
    }
    else {
      
      length.exon=exon.cor[mut.exon[[1]]]-as.numeric(variants.features.fr$coding.pos[x])
      
      
    }
  }
  
})

variants.features.fr$PTC.2.EJC <- lapply(variants.features.fr$PTC.2.EJC, function(x) if (is.null(x)) NA else x)
variants.features.fr$PTC.2.EJC <- unlist(variants.features.fr$PTC.2.EJC)

### PTC.2.end
variants.features.fr$PTC.2.end <-  sapply(1:nrow(variants.features.fr),function(x)
  
{
  variants.features.fr$cds_length[x]-as.numeric(variants.features.fr$coding.pos[x])
  
})

### PTC.2.start

variants.features.fr$PTC.2.start <-  sapply(1:nrow(variants.features.fr),function(x)
  
{
  as.numeric(variants.features.fr$coding.pos[x])
  
})

### last exon rule

variants.features.fr$last.exon <-  sapply(1:nrow(variants.features.fr),function(x)
  
{
  print(x)
  if(!is.na(variants.features.fr$cds_exons[x])& !variants.features.fr$mut.exon[x]=='NA' & !is.null(variants.features.fr$mut.exon[x][[1]])){
    if(variants.features.fr$exon_count[x]=='1') {
      
      paste('NA')
    } else if (variants.features.fr$exon_count[x]==variants.features.fr$mut.exon[x]) {
      
      paste('lastexon')
    }
    else {
      
      paste('notlastexon')
    }
    
    
    
  }
  
  
})

### penultimate exon rule

variants.features.fr$penultimate.exon <-  sapply(1:nrow(variants.features.fr),function(x)
  
{
  
  print(x)
  if(!is.na(variants.features.fr$cds_exons[x]) & !variants.features.fr$mut.exon[x]=='NA' & !is.null(variants.features.fr$mut.exon[x][[1]])) {
    if(as.numeric(variants.features.fr$exon_count[x])==1) {
      
      paste('not penultimate.last50bp')
    } else if ( as.numeric(variants.features.fr$exon_count[x])>1 & as.numeric(variants.features.fr$exon_count[x])==variants.features.fr$mut.exon[x][[1]]+1) {
      
      exon.cor <- as.numeric(unlist(strsplit(as.character(variants.features.fr$cds_exons[x]),',')))
      if(length (exon.cor)>=2){
        fifty.pos <- exon.cor[length(exon.cor)-1]-50
        penultimate.exon.length <- exon.cor[length(exon.cor)-1]-exon.cor[length(exon.cor)-2]
      }
      if(length (exon.cor)==2){
        fifty.pos <- exon.cor[length(exon.cor)-1]-50
        penultimate.exon.length <- exon.cor[1]
      }
      if (as.numeric(variants.features.fr$coding.pos[x])>=fifty.pos & as.numeric(variants.features.fr$coding.pos[x])<=exon.cor[length(exon.cor)-1] & penultimate.exon.length>=50 ){
        paste('penultimate.last50bp')
        
      } else {
        paste('not penultimate.last50bp')
        
      }
      
    } else {
      
      paste('not penultimate.last50bp')
      
    }
    
  } else {
    
    paste('not penultimate.last50bp')
  }
  
})


######## get the unique variants look at the distinguishing features####

variants.features.fr$ALLELE.RAT<- 1-variants.features.fr$altCount/variants.features.fr$totalCount
variants.features.fr$NMD.ESCAPEE <- rep('NA',nrow(variants.features.fr))
esc.ind <- which(variants.features.fr$ALLELE.RAT<=0.65 & variants.features.fr$ALLELE.RAT>=0.35)
nonesc.ind <-  which(variants.features.fr$ALLELE.RAT>0.65)
variants.features.fr[esc.ind,'NMD.ESCAPEE'] <- 'TRUE'
variants.features.fr[nonesc.ind,'NMD.ESCAPEE'] <- 'FALSE'

#### get the distinct variants 
### not to count any variant more than once 

tissue.var <- variants.features.fr
uni.var <- unique(tissue.var$key)

freq.fr <- data.frame(table(tissue.var$key))

sub.ind <- NULL
for (x in 1:length(uni.var))
{
  
  var <- uni.var[x]
  ind <- which(tissue.var$key==var)
  ind.1 <- sample(length(ind),1)
  ind.2 <- ind[ind.1]
  sub <- tissue.var[ind,]

  
  sub.ind <- c(sub.ind,ind.2)
}

tissue.var.sub <- tissue.var[sub.ind,]

df <- tissue.var.sub

freq.fr$key <- freq.fr$Var1

df <- merge(df,freq.fr,by='key')
#only stopgain variants
#df.sub <- df[which(df$V2=='stopgain'),]

#checking if multi-alleles
df.sub$REF_len <- nchar(df.sub$REF_ALLELE)
df.sub$ALT_len <- nchar(df.sub$ALT_ALLELE)

head(df)
df_snps <- df[df$REF_len == 1 & df$ALT_len == 1, ]


                                         
