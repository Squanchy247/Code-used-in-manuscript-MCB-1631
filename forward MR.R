install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
library(ggplot2)
library(TwoSampleMR)
library(foreach)

setwd("")
library(data.table)
MG<-fread(input = "finngen_R11_K11_ALCOLIV.gz",sep = "\t",header = TRUE)#!!!!!!!!!!!!!!!!!!!!!!!!!!
head(MG)
dbm=read.table("cathepsin_id.txt",header = T,sep = "\t")

sbdbm=as.vector(dbm$id)
foreach(i=sbdbm, .errorhandling = "pass") %do%{
  expo_rt<-read.table(file =paste0(i,".clump.txt"),header = T,sep = "\t")

  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = "finngen_R11_K11_ALCOLIV.gz",#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sep = "\t",
    snp_col = "rsids",
    phenotype_col = "phenotype",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval")

harmonise_rt <- harmonise_data(
  exposure_dat =  expo_rt, 
  outcome_dat = outc_rt,action=2)

harmonise_rt$R2 <- (2 * (harmonise_rt$beta.exposure^2) * harmonise_rt$eaf.exposure * (1 - harmonise_rt$eaf.exposure)) /
  (2 * (harmonise_rt$beta.exposure^2) * harmonise_rt$eaf.exposure * (1 - harmonise_rt$eaf.exposure) +
     2 * harmonise_rt$samplesize.exposure*harmonise_rt$eaf.exposure * (1 - harmonise_rt$eaf.exposure) * harmonise_rt$se.exposure^2)
harmonise_rt$f <- harmonise_rt$R2 * (harmonise_rt$samplesize.exposure - 2) / (1 - harmonise_rt$R2)
harmonise_rt$meanf<- mean( harmonise_rt$f)
harmonise_rt<-harmonise_rt[harmonise_rt$f>10,]

forward_result<- mr(harmonise_rt)
result_or=generate_odds_ratios(forward_result) 
dir.create(i) 
filename=i
write.table(harmonise_rt, file =paste0(filename,"/harmonise.txt"),row.names = F,sep = "\t",quote = F)
write.table(result_or[,5:ncol(result_or)],file =paste0(filename,"/OR.txt"),row.names = F,sep = "\t",quote = F)
pleiotropy=mr_pleiotropy_test(harmonise_rt)
write.table(pleiotropy,file = paste0(filename,"/pleiotropy.txt"),sep = "\t",quote = F)
heterogeneity=mr_heterogeneity(harmonise_rt)
write.table(heterogeneity,file = paste0(filename,"/heterogeneity.txt"),sep = "\t",quote = F)
presso=run_mr_presso(harmonise_rt,NbDistribution = 1000)
capture.output(presso,file = paste0(filename,"/presso.txt"))

singlesnp_res<- mr_singlesnp(harmonise_rt)
singlesnpOR=generate_odds_ratios(singlesnp_res)
write.table(singlesnpOR,file=paste0(filename,"/singlesnpOR.txt"),row.names = F,sep = "\t",quote = F)
}
