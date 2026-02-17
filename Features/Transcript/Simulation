df <- df[df$Freq.cat=='Rare/Common variants',]

df.sub.new <- df
df.sub.new.1 <- NULL

#df.sub.new.1 <- NULL
### 100 simulations

for (t in 1:100) {
  
  print(t)
  df.sub.new <- df
  for ( i in 1:nrow(df.sub.new)){
    print(i)
    ind <- which(fr.var.can$variantID%in%df.sub.new$variantID[i])
    ind.1 <- sample(ind,1, replace=F) 
    new.refcount <- fr.var.can$refCount[ind.1]
    new.totalcount <- fr.var.can$totalCount[ind.1]
    new.ALLELE.RAT <- new.refcount / new.totalcount
    df.sub.new$ALLELE.RAT[i] <-  new.ALLELE.RAT
    
  }
  #df.sub.new.changed <- df.sub.new
  df.sub.new.1 <- rbind(df.sub.new,df.sub.new.1)
  
}
df.sub <- df.sub.new.1 %>% group_by(key) %>% mutate(ALLELE.RAT = median(ALLELE.RAT)) %>% unique()

hist(df.sub$ALLELE.RAT, main = "Rare/Common (after simulation)")
#optional
#hist(df.sub$ALLELE.RAT, main = "Rare/Common (before simulation)")

#merge the Rare/common variants back to the main data frame:
df$ALLELE.RAT[df$Freq.cat == "Rare/Common variants"] <- df.sub$ALLELE.RAT
summary(df$ALLELE.RAT[df$Freq.cat == "Rare/Common variants"])
