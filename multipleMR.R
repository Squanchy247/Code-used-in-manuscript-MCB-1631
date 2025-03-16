install.packages("MendelianRandomization")

library(MendelianRandomization)
library(TwoSampleMR)

setwd("")

df=read.table("cathepsin_id.txt",header = T,sep = "\t")
id=as.vector(df$id)
exposure_id=c(id)
exposure_rt<- mv_extract_exposures(exposure_id,pval_threshold = 5e-06)
outcome_dat <- read
outcome_dat <- read_outcome_data(
  snps = exposure_rt$SNP,
  filename = "finngen_R11_K11_APPENDACUT.gz",
  sep = "\t",
  snp_col = "rsids",
  phenotype_col = "phenotype",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval")

mvhr_rt <- mv_harmonise_data(exposure_rt, outcome_dat)


mvmrdata <- mr_mvinput(bx = mvhr_rt$exposure_beta, bxse = mvhr_rt$exposure_se, by = mvhr_rt$outcome_beta, byse = mvhr_rt$outcome_se, correlation =matrix())


mvmrdata <- mr_mvinput(bx = mvhr_rt$exposure_beta, bxse = mvhr_rt$exposure_se, by = mvhr_rt$outcome_beta, byse = mvhr_rt$outcome_se, correlation =matrix())

#res <- mv_multiple(mvhr_rt)
#____________________________________________________________________________
res_ivw<-mr_mvivw(mvmrdata)

exposure_names <- names(res_ivw@Estimate)
estimates <- res_ivw@Estimate
std_errors <- res_ivw@StdError
ci_lower <- res_ivw@CILower
ci_upper <- res_ivw@CIUpper
p_values <- res_ivw@Pvalue

ors <- exp(estimates)
lci_95_ors <- exp(ci_lower)
uci_95_ors <- exp(ci_upper)

results_df <- data.frame(
  id.exposure = paste("id:", exposure_names, sep = ""),
  exposure = exposure_names,
  id.outcome = "id:outcome",
  outcome = "outcome",
  nsnp = 101,
  b = estimates,
  se = std_errors,
  pval = p_values,
  or = ors,
  lci95 = lci_95_ors,
  uci95 = uci_95_ors,
  stringsAsFactors = FALSE
)

write.table(results_df, file = "res_ivw_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#——————————————————————————————————————————————————————————————————
res_fe_ivw <- mr_mvivw(mvmrdata, model = "fixed")

exposure_names <- names(res_fe_ivw@Estimate)
estimates <- res_fe_ivw@Estimate
std_errors <- res_fe_ivw@StdError
ci_lower <- res_fe_ivw@CILower
ci_upper <- res_fe_ivw@CIUpper
p_values <- res_fe_ivw@Pvalue

ors <- exp(estimates)
lci_95_ors <- exp(ci_lower)
uci_95_ors <- exp(ci_upper)

results_df <- data.frame(
  id.exposure = paste("id:", exposure_names, sep = ""),
  exposure = exposure_names,
  id.outcome = "LVt7sv",
  outcome = "outcome",
  nsnp = 1,
  b = estimates,
  se = std_errors,
  pval = p_values,
  or = ors,
  lci95 = lci_95_ors,
  uci95 = uci_95_ors,
  stringsAsFactors = FALSE
)

write.table(results_df, file = "res_fe_ivw_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#——————————————————————————————————————————————————————————————————
res_re_ivw <- mr_mvivw(mvmrdata, model = "random")

exposure_names <- names(res_re_ivw@Estimate)
estimates <- res_re_ivw@Estimate
std_errors <- res_re_ivw@StdError
ci_lower <- res_re_ivw@CILower
ci_upper <- res_re_ivw@CIUpper
p_values <- res_re_ivw@Pvalue

id_exposure <- paste("id:", exposure_names, sep="")
id_outcome <- "LVt7sv" 

results_df <- data.frame(
  id.exposure = id_exposure,
  exposure = exposure_names,
  id.outcome = id_outcome,
  outcome = "outcome",
  nsnp = 101,
  b = estimates,
  se = std_errors,
  pval = p_values,
  or = exp(estimates),
  lci95 = exp(ci_lower),
  uci95 = exp(ci_upper),
  stringsAsFactors = FALSE
)

write.table(results_df, file = "res_re_ivw_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#MR-Egger
#______________________________________________________________________________
res_egger<-mr_mvegger(mvmrdata)

exposure_names <- res_egger@Exposure
estimates <- res_egger@Estimate
std_errors <- res_egger@StdError.Est
ci_lower <- res_egger@CILower.Est
ci_upper <- res_egger@CIUpper.Est
p_values <- res_egger@Pvalue.Est

ors <- exp(estimates)
lci_95_ors <- exp(ci_lower)
uci_95_ors <- exp(ci_upper)

results_df <- data.frame(
  id.exposure = paste("id:", exposure_names, sep = ""),
  exposure = exposure_names,
  id.outcome = "id:outcome",
  outcome = "outcome",
  nsnp = 101,
  b = estimates,
  se = std_errors,
  pval = p_values,
  or = ors,
  lci95 = lci_95_ors,
  uci95 = uci_95_ors,
  stringsAsFactors = FALSE
)

write.table(results_df, file = "res_egger_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")
#Median
#______________________________________________________________________________
res_median<-mr_mvmedian(mvmrdata)

exposure_names <- res_median@Exposure
estimates <- res_median@Estimate
std_errors <- res_median@StdError
ci_lower <- res_median@CILower
ci_upper <- res_median@CIUpper
p_values <- res_median@Pvalue

ors <- exp(estimates)
lci_95_ors <- exp(ci_lower)
uci_95_ors <- exp(ci_upper)

results_df <- data.frame(
  id.exposure = paste("id:", exposure_names, sep = ""),
  exposure = exposure_names,
  id.outcome = "id:outcome",
  outcome = "outcome",
  nsnp = 101,
  b = estimates,
  se = std_errors,
  pval = p_values,
  or = ors,
  lci95 = lci_95_ors,
  uci95 = uci_95_ors,
  stringsAsFactors = FALSE
)

write.table(results_df, file = "res_median_results.txt", row.names = FALSE, quote = FALSE, sep = "\t")