library(TwoSampleMR)
library(ggplot2)
library(foreach)

setwd("")

IDgo=read.table("cathepsin_id.txt",header = T,sep = "\t")

eachdbm=as.vector(IDgo$id) 

foreach(i=eachdbm, .errorhandling = "pass") %do%{
  
  expo_rt<-read.table("expo_rt_out_APPENDACUT.txt",header = T,sep = "\t")#！ Don't forget to use the file generated in step1 named "expo_rt_out_XXX"
  expo_rt$samplesize.exposure=451025
  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = paste0(i,".txt.gz"),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "eaf",
    pval_col = "p"
  )
  
  
  harmonise_rt <- harmonise_data(
    exposure_dat =  expo_rt, 
    outcome_dat = outc_rt,action=2)
  mr_result<- mr(harmonise_rt)
  
  harmonise_rt$R2 <- (2 * (harmonise_rt$beta.exposure^2) * harmonise_rt$eaf.exposure * (1 - harmonise_rt$eaf.exposure)) /
    (2 * (harmonise_rt$beta.exposure^2) * harmonise_rt$eaf.exposure * (1 - harmonise_rt$eaf.exposure) +
       2 * harmonise_rt$samplesize.exposure*harmonise_rt$eaf.exposure * (1 - harmonise_rt$eaf.exposure) * harmonise_rt$se.exposure^2)
  harmonise_rt$f <- harmonise_rt$R2 * (harmonise_rt$samplesize.exposure - 2) / (1 - harmonise_rt$R2)
  harmonise_rt$meanf<- mean(harmonise_rt$f)
  harmonise_rt<-harmonise_rt[harmonise_rt$f>10,]
  
  
  result_or=generate_odds_ratios(mr_result) 
  
  if(result_or$pval[3]<1){
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
    
  }
}

