library(data.table)
#load the dataset
variant_anno<- fread('~/TOPMed_v1_stopgain_frameshift_gencode_v38.exonic_variant_function')
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

variant_anno <- variant_anno[c(frame.ind,stop.ind),]

variant_anno.fr <- variant_anno[frame.ind,]

for(i in 1:nrow(variant_anno.fr))
{
    print(i)
    
      
      if(length(which(fr.var.can$position==variant_anno.fr$V5[i]-1))>0 | length(which(fr.var.can$position==variant_anno.fr$V5[i]))>0 ) {
          
           ind <- which(fr.var.can$position==variant_anno.fr$V5[i]-1 & fr.var.can$contig==variant_anno.fr$V4[i])
           ind.1 <- which(fr.var.can$position==variant_anno.fr$V5[i] & fr.var.can$contig==variant_anno.fr$V4[i])
           ind.2 <- c(ind,ind.1)
           for ( x in ind.2) 
           {
              if(nchar(fr.var.can$REF_ALLELE[x])> nchar(fr.var.can$ALT_ALLELE[x])){
                #chr <- mapply(setdiff, strsplit(fr.var.can$REF_ALLELE[x], "\\s+"), strsplit(fr.var.can$ALT_ALLELE[x], "\\s+"))
                chr <- substr(fr.var.can$REF_ALLELE[x],nchar(fr.var.can$ALT_ALLELE[x])+1,nchar(fr.var.can$REF_ALLELE[x]))
                fr.var.can$key[x]= paste(fr.var.can$contig[x],':',fr.var.can$position[x]+1,'_',chr,'>','-',sep='')
                
              } else if (nchar(fr.var.can$REF_ALLELE[x])< nchar(fr.var.can$ALT_ALLELE[x])) {
                #chr <- mapply(setdiff, strsplit(fr.var.can$ALT_ALLELE[x], "\\s+"), strsplit(fr.var.can$REF_ALLELE[x], "\\s+"))
                chr <- substr(fr.var.can$ALT_ALLELE[x],nchar(fr.var.can$REF_ALLELE[x])+1,nchar(fr.var.can$ALT_ALLELE[x]))
                fr.var.can$key[x]= paste(fr.var.can$contig[x],':',fr.var.can$position[x],'_','-','>',chr,sep='')
                
              }
                
            }
             
             
             
          }
        
        
}
      
    
  
  
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
