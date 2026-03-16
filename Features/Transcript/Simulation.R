#df <- df[df$Freq.cat=='Rare/Common variants',]

df.sub.new <- df
df.sub.new.1 <- NULL

#df.sub.new.1 <- NULL
### 100 simulations
missing_variants <- 0

for (t in 1:100) {
  
  print(t)
  df.sub.new <- df
  for ( i in 1:nrow(df.sub.new)){
    print(i)
    #ind <- which(fr.var.can$variantID%in%df.sub.new$variantID[i]) #fr.var.can file is made at the beginning check Data tab
    ind <- which(fr.var.can$variantID == df.sub.new$variantID[i])
    
    if(length(ind)==0){
      print(df.sub.new$variantID[i])
      missing_variants <- missing_variants + 1
      df.sub.new$ALLELE.RAT[i] <- NA

      next
    }
  
    
    ind.1 <- sample(ind,1) 
    new.refcount <- fr.var.can$refCount[ind.1]
    new.totalcount <- fr.var.can$totalCount[ind.1]
    new.ALLELE.RAT <- new.refcount / new.totalcount
    df.sub.new$ALLELE.RAT[i] <-  new.ALLELE.RAT
    
  }
  #df.sub.new.changed <- df.sub.new
  df.sub.new.1 <- rbind(df.sub.new,df.sub.new.1)
  
}
df.sub <- df.sub.new.1 %>% group_by(key) %>% mutate(ALLELE.RAT = median(ALLELE.RAT, na.rm=TRUE)) %>% unique()

hist(df.sub$ALLELE.RAT, main = "Rare/Common (after simulation)")
#optional
#hist(df.sub$ALLELE.RAT, main = "Rare/Common (before simulation)")

#merge the Rare/common variants back to the main data frame:
df$ALLELE.RAT[df$Freq.cat == "Rare/Common variants"] <- df.sub$ALLELE.RAT
summary(df$ALLELE.RAT[df$Freq.cat == "Rare/Common variants"])
