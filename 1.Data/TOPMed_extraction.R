#libraries
library(data.table)
library(ggplot2)
library(gridExtra)
library(plyr)
library(dplyr)

# ASE EXTRACTION PIPELINE FOR PROTEIN-TRUNCATING VARIANTS (PTVs)
# TOPMed Freeze 9b
#
# Input:
#   - TOPMed Freeze 9b VCF
#   - matched RNA-seq BAM files
#
# Output:
#   - deduplicated PTV VCF
#   - ASEReadCounter allele counts

bcftools view \
    --drop-genotypes \
    -i 'INFO/ANN ~ "frameshift_variant" || INFO/ANN ~ "stop_gained"' \
    TOPMed_Freeze9b.vcf.gz \
    -Oz \
    -o TOPMed_PTVs.vcf.gz

#

bcftools norm \
    --rm-dup all \
    TOPMed_PTVs.vcf.gz \
    -Oz \
    -o TOPMed_PTVs.dedup.vcf.gz
#
tabix -p vcf TOPMed_PTVs.dedup.vcf.gz

gatk ASEReadCounter \
    -R GRCh38.fa \
    -I sample.rna.bam \
    -V TOPMed_PTVs.dedup.vcf.gz \
    -O sample_ASE_counts.table \
    --min-depth 1 \
    --min-mapping-quality 255 \
    --min-base-quality 10
###############################################################################
# Summary
#
# Initial extracted variants:
#   77,918 stop_gained + frameshift variants
#
# After duplicate removal:
#   77,408 variants
#
# ASE counts computed across:
#   11,630 TOPMed donors with matched RNA-seq
###############################################################################
#Genotype and ASE files were merged together
setwd('~/ASE_genotype')

files <- list.files('.','.vcf')
het_get <- data.frame()
for (i in 1: length(files)){
print(i)
TOR <- read.table(files[i],header=TRUE)
ind <- which(TOR[,ncol(TOR)]=='0/1')
TOR.sub <- TOR[ind,]
colnames(TOR.sub)[ncol(TOR.sub)]<-'sample'
het_get <- rbind(het_get, TOR.sub)
}
# Combine the individual data frames
het_get <- do.call(rbind, het_get_list)

write.table(het_get, file='~/merged_het_get.txt', row.names=FALSE, sep="\t", quote=FALSE)

fr.var.can <- het_get

fr.var.can$CHROM <- het_get$contig
fr.var.can$POS <- het_get$position
fr.var.can$REF_ALLELE <- sapply(1:nrow(fr.var.can),function(x){
  strsplit(fr.var.can$variantID[x],'_')[[1]][3]
})

fr.var.can$ALT_ALLELE <- sapply(1:nrow(fr.var.can),function(x){
  
  strsplit(fr.var.can$variantID[x],'_')[[1]][4]
  
})
fr.var.can$key <- paste(fr.var.can$CHR,':',fr.var.can$POS,'_',fr.var.can$REF_ALLELE,'>',fr.var.can$ALT_ALLELE,sep='')


#### annotate TOPMed  data with annovar

fr.var.can$ZYG <- rep('Het',nrow(fr.var.can))
fr.var.can$FILTER <- rep('PASS',nrow(fr.var.can))
ID <- rep('.',nrow(fr.var.can))
merged <- cbind(fr.var.can$CHR,fr.var.can$POS,fr.var.can$POS, as.character(fr.var.can$REF_ALLELE),as.character(fr.var.can$ALT_ALLELE),ID, as.character(fr.var.can$ZYG),as.character(fr.var.can$FILTER))

write.table(merged,'~/NMD_TOPMed/TOPMed_v1_stopgain_frameshift.vcf',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
system(paste('perl ~/annovar/convert2annovar.pl -format vcf4 ~/NMD_TOPMed/TOPMed_v1_stopgain_frameshift.vcf > ~/NMD_TOPMed/TOPMed_v1_stopgain_frameshift.avinput',sep=''))

system(paste('perl ~/annovar/annotate_variation.pl -build hg38 -out ~/NMD_TOPMed/TOPMed_v1_stopgain_frameshift_gencode_v38 -dbtype ensGene ~/NMD_TOPMed/TOPMed_v1_stopgain_frameshift.avinput ~/annovar/tempdir',sep=''))

