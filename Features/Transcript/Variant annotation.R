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
variants.features.fr$PTC.2.EJC <- sapply(1:nrow(variants.features.fr), function(x) {

  if (!is.na(variants.features.fr$cds_exons[x]) &&
      !is.na(variants.features.fr$mut.exon[x])) {

    exon.cor <- as.numeric(unlist(strsplit(
      as.character(variants.features.fr$cds_exons[x]), ',')))

    mut.exon <- variants.features.fr$mut.exon[x]

    return(exon.cor[mut.exon] - as.numeric(variants.features.fr$coding.pos[x]))
  }

  return(NA_real_)
})


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
variants.features.fr$last.exon <- sapply(1:nrow(variants.features.fr), function(x) {

  if (!is.na(variants.features.fr$cds_exons[x]) &&
      !is.na(variants.features.fr$mut.exon[x])) {

    if (as.numeric(variants.features.fr$exon_count[x]) == 1) {
      return(NA_character_)
    } else if (as.numeric(variants.features.fr$exon_count[x]) ==
               variants.features.fr$mut.exon[x]) {
      return("lastexon")
    } else {
      return("notlastexon")
    }
  }

  return(NA_character_)
})      

### penultimate exon rule

variants.features.fr$penultimate.exon <- sapply(1:nrow(variants.features.fr), function(x) {

  if (!is.na(variants.features.fr$cds_exons[x]) &&
      !is.na(variants.features.fr$mut.exon[x])) {

    exon_count <- as.numeric(variants.features.fr$exon_count[x])
    mut_exon <- variants.features.fr$mut.exon[x]

    # single exon → not applicable
    if (exon_count == 1) {
      return("not penultimate.last50bp")
    }

    # check penultimate exon
    if (exon_count > 1 && exon_count == mut_exon + 1) {

      exon.cor <- as.numeric(unlist(strsplit(
        as.character(variants.features.fr$cds_exons[x]), ',')))

      if (length(exon.cor) >= 2) {

        penultimate_end <- exon.cor[length(exon.cor) - 1]
        penultimate_start <- ifelse(length(exon.cor) > 2,
                                   exon.cor[length(exon.cor) - 2],
                                   0)

        penultimate_length <- penultimate_end - penultimate_start
        fifty_pos <- penultimate_end - 50

        coding_pos <- as.numeric(variants.features.fr$coding.pos[x])

        if (coding_pos >= fifty_pos &&
            coding_pos <= penultimate_end &&
            penultimate_length >= 50) {

          return("penultimate.last50bp")
        }
      }
    }

    return("not penultimate.last50bp")
  }

  return(NA_character_)
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


                                         
