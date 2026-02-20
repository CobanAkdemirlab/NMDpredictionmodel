
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


