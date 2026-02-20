df.sub <- df[which(df$V2=='stopgain'),]

##### get the PTC codon identity for only stopgain variants

df.sub$PTC.id <- NULL
for ( i in 1:nrow(df.sub)){
  print(i)
  ind <- which(names(cds_seqs)==df.sub$txnames[i])
  ens.annotated <- df.sub$V3[i]
  pos <- strsplit(ens.annotated,':|,')
  cds.pos <- pos[[1]][grep(df.sub$txnames[i],pos[[1]])+2]
  changed.aa <- substr(cds.pos,nchar(cds.pos),nchar(cds.pos))
  
  if(length(grep('c.',cds.pos))==1) {
    cds.pos <- as.numeric(gsub("\\D", "", cds.pos))
    if(cds.pos%%3==1){
      
      PTC.id <- substr(as.character(cds_seqs[ind][1]),cds.pos,cds.pos+2)
      substr(PTC.id,1,1) <- changed.aa
      df.sub$PTC.id[i] <- PTC.id
      
    } else if(cds.pos%%3==2){
      
      PTC.id <- substr(as.character(cds_seqs[ind][1]),cds.pos-1,cds.pos+1)
      substr(PTC.id,2,2) <- changed.aa
      df.sub$PTC.id[i] <- PTC.id
      
    } else {
      PTC.id <- substr(as.character(cds_seqs[ind][1]),cds.pos-2,cds.pos)
      substr(PTC.id,3,3) <- changed.aa
      df.sub$PTC.id[i] <- PTC.id
      
    }
    
  } else 
  {
    'NA'
    
  }
  
  
  
}
df.sub$PTC.id <- as.factor(df.sub$PTC.id)
#### ### get the rare/common variants

df.sub <- df.sub[which(df.sub$Freq>1),]

