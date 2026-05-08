##load('TOPMed_whlbld_variantsusedfortheanalysis.RData')

# Ensure df$mut.exon is numeric
df$mut.exon <- as.numeric(df$mut.exon)

# Calculate EJC.downstream after making sure mut.exon is numeric
df$EJC.downstream <- df$exon_count - df$mut.exon
#categorization
df$EJC.downstream.cut <- cut(
    df$EJC.downstream,
    breaks = c(-Inf, 2, 5, 9, 14, Inf),
    labels = c("0-2", "3-5", "6-9", "10-14", ">14"),
    right = TRUE
)

table(df$EJC.downstream.cut)

#### exception 1 ### the ones in the last exon or penutimate.last50bp


exc.1 <- which(df$last.exon=='lastexon')
exc.2 <- which(df$penultimate.exon=='penultimate.last50bp')


#### exception 1 ### merge the ones in the last exon and penultimate exon


exc.3 <- which(df$last.exon=='lastexon' | df$penultimate.exon=='penultimate.last50bp')







#### make the canonical rule feature

df$last.EJC <- rep('NA',nrow(df))
df$last.EJC[exc.1] <- 'last.exon'
df$last.EJC[exc.2] <- 'penultimate.last50bp'
df$last.EJC[-c(exc.1,exc.2)] <- 'upstream'

df$downstream <- rep('NA',nrow(df))

df$downstream[exc.3] <- 'downstream of last EJC'
df$downstream[-exc.3] <- 'upstream of last EJC'


### make the Freq cat

df$Freq.cat <- rep('NA',nrow(df))
ind <- which(df$Freq==1)
df$Freq.cat[ind] <- 'Ultra-rare variants'
ind <- which(df$Freq>1)
df$Freq.cat[ind] <- 'Rare/Common variants'

table(df$Freq.cat)
