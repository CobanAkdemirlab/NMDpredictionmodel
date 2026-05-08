library(data.table)
#load the dataset (only snv and is ptc)
variant_anno<- fread('~/NMD-TOPMed/TOPMed_v1_stopgain_gencode_v38.exonic_variant_function')
### start of annotation
variant_anno$key <- paste(variant_anno$V4,':',variant_anno$V5,'_',variant_anno$V7,'>',variant_anno$V8,sep='')

#### merge variant annotations with the genotype data frame
variant_anno <- data.frame(variant_anno)



#### Remove the nonunique columns
variant_anno <- variant_anno[,c(2:8,12)]
variant_anno <- unique(variant_anno)

### Upload the Variant Annotation Data
load('~/fr.var.can.RData')

## get only the stopgains and frameshifts and then get the key from fr.var.can 
frame.ind <- grep('frameshift',variant_anno$V2)

stop.ind <- grep ('stopgain',variant_anno$V2)

variant_anno <- variant_anno[c(stop.ind),]

      
    
  
  
# Merge both datasets
variant_anno_merged <- merge(variant_anno,fr.var.can,by='key')


#### annotate them with NMD features
variant.anno.hiqual <- variant_anno_merged


library(biomaRt)
mart <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
BM.info <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol","transcript_is_canonical"),mart=mart)
BM.info.can <- BM.info[which(BM.info$transcript_is_canonical=='1'),]

### Get the canonical transcripts from BioMart
txnames <-  sapply(1:nrow(variant.anno.hiqual),function(x)
  
{
  print(x)
  ens.annotated <- variant.anno.hiqual$V3[x]
  gene <- strsplit(ens.annotated,':')[[1]][1]
  can.transcript <- BM.info.can[which(BM.info.can$hgnc_symbol==gene),'ensembl_transcript_id']
  if(length(can.transcript)>0){
  all.tx <- strsplit(ens.annotated,':')
  if(length(grep(can.transcript,all.tx[[1]]))==1){
  all.tx[[1]][grep(can.transcript,all.tx[[1]])]
  } else{
    
    'NA'
  }
  } else {
    
    'NA'
  }
  
  
})
variant.anno.hiqual$txnames <- txnames
### only get the annotations with canonical transcripts
variant.anno.hiqual <- variant.anno.hiqual[which(variant.anno.hiqual$txnames!='NA'),]
#save
#add gencode V26

