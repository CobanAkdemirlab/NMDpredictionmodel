#####
library(parallel)
library(dplyr)
df.sub.new <- df
df.sub.new.1 <- NULL

#df.sub.new.1 <- NULL
### 100 simulations
missing_variants <- 0
ind <- which(fr.var.can$totalCount<8)
fr.var.can.new <- fr.var.can[-ind,]
load('~/Downloads/TOPMed_stopgain_fulldf_feb1826.RData')

### these ar 6,019 variants
common.vrnts <- intersect(df$key,fr.var.can.new$key)
## subset of fr.var.can.new and df on the common vrnts
fr.var.can.new <- fr.var.can.new[which(fr.var.can.new$key%in%common.vrnts),]
df <- df[which(df$key%in%common.vrnts),]

### 


fr <- data.frame(table(fr.var.can.new$variantID))
ultra <- unique(fr[which(fr$Freq==1),'Var1'])
common <- unique(fr[which(fr$Freq>1),'Var1'])

df.rare <- df[which(df$variantID%in%ultra),]
df.common <- df[which(df$variantID%in%common),]


### first

for ( i in 1:length(ultra))
{
  ind <- which(fr.var.can.new$variantID%in%ultra[i])
  df.rare[which(df.rare$variantID%in%ultra[i]),'refCount'] <- fr.var.can.new[ind,'refCount']
  df.rare[which(df.rare$variantID%in%ultra[i]),'altCount'] <- fr.var.can.new[ind,'altCount']
  df.rare[which(df.rare$variantID%in%ultra[i]),'totalCount'] <- fr.var.can.new[ind,'totalCount']
  df.rare[which(df.rare$variantID%in%ultra[i]),'ALLELE.RAT'] <- df.rare[which(df.rare$variantID%in%ultra[i]),'refCount']/df.rare[which(df.rare$variantID%in%ultra[i]),'totalCount'] 
  
}

### 

df.sub.new <- df.common
df.sub.new.1 <- NULL



for (t in 1:100) {
  
  print(t)
  df.sub.new <- df.common
  for ( i in 1:nrow(df.sub.new)){
    print(i)
    #ind <- which(fr.var.can$variantID%in%df.sub.new$variantID[i]) #fr.var.can file is made at the beginning check Data tab
    ind <- which(fr.var.can.new$variantID == df.sub.new$variantID[i])
    
  
    
    
    ind.1 <- sample(ind,1) 
    new.refcount <- fr.var.can.new$refCount[ind.1]
    new.totalcount <- fr.var.can.new$totalCount[ind.1]
    new.ALLELE.RAT <- new.refcount / new.totalcount
    df.sub.new$ALLELE.RAT[i] <-  new.ALLELE.RAT
    
  }
  #df.sub.new.changed <- df.sub.new
  df.sub.new.1 <- rbind(df.sub.new,df.sub.new.1)
  
}

library(dplyr)

## new simulated data frame for common variants
df.sim.common <- df.sub.new.1 %>% group_by(key) %>% mutate(ALLELE.RAT = median(ALLELE.RAT, na.rm=TRUE)) %>% distinct(.keep_all = TRUE)



## merge common and rare data and this is the new data frame

df.sim <- rbind(df.rare,df.sim.common)

### please define NMD.escape since allele.rat change

df.sim$NMD.ESCAPEE <- rep('NA',nrow(df.sim))
esc.ind <- which(df.sim$ALLELE.RAT<=0.65 & df.sim$ALLELE.RAT>=0.35)
nonesc.ind <-  which(df.sim$ALLELE.RAT>0.65)
df.sim[esc.ind,'NMD.ESCAPEE'] <- 'TRUE'
df.sim[nonesc.ind,'NMD.ESCAPEE'] <- 'FALSE'
