#!/usr/bin/Rscript
  
af <- read.table("af_101SNPs_80inds.txt", sep = "\t")

af_mod <- af
af_mod[af>0.7] <- 0
af_mod[af>0.5 & af<0.7] <- 0.5


write.table(af_mod, "afmod_101SNPs_80inds.txt", sep = "\t", col.names = F, row.names = F)
