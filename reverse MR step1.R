setwd("")

library(data.table)
source("ld_clump.R")
source("ld_matrix.R")
source("backwards.R")

expo_rt=fread("finngen_R11_CHRONSMALL.gz", header = T)


expo_rt=expo_rt[expo_rt$pval<5e-6,]

expo_rt2=expo_rt[,c("rsids","pval")]
colnames(expo_rt2)=c("rsid", "pval")


clumdf <- ld_clump_local(dat = expo_rt2, clump_kb = 10000, clump_r2 = 0.001,clump_p=1,
                     bfile ="D:\\data_maf0.01_rs_ref", 
                     plink_bin = "D:\\plink.exe")

expo_rt3=expo_rt[which(expo_rt$rsids%in%clumdf$rsid),]

write.table(expo_rt3,"expo_rt_CHRONSMALL111.txt",row.names = F,sep = "\t",quote = F)

####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!wait till expo_rt_CHRONSMALL111.txt comes out
library(TwoSampleMR)

df=read_exposure_data(filename = "expo_rt_CHRONSMALL111.txt"
                      sep = "\t",
                      snp_col = "rsids",
                      phenotype_col = "phenotype",
                      beta_col = "beta",
                      se_col = "sebeta",
                      effect_allele_col = "alt",
                      other_allele_col = "ref",
                      eaf_col = "af_alt",
                      pval_col = "pval")

write.table(df,"expo_rt_out_CHRONSMALL1111.txt",row.names = F,sep = "\t",quote = F)
# then copy expo_rt_out_CHRONSMALL1111.txt to another fileclip:step2


