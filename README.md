# Code-used-in-manuscript-MCB-1631
Code used in manuscript “Genetic and biomechanics insights into cathepsins and non-cancerous digestive diseases: a bidirectional two-sample Mendelian randomization study”

1.Forward MR (forward MR.R)
Set working directory using setwd("").
Read outcome data from the "outcome file" and exposure IDs from cathepsin_id.txt.
Loop through each exposure ID: read exposure data, harmonize with outcome data, calculate statistics (e.g., R2, f - statistic), conduct MR analysis, and save results in an ID - named directory.


2.Reverse MR

2.1 Step 1 (reverse MR step1.R)
Set working directory with setwd("").
Source ld_clump.R, ld_matrix.R, and backwards.R.
Read and filter initial exposure data from the outcome file (p - value < 5e - 6).
Perform LD clumping with PLINK (clump_kb = 10000, clump_r2 = 0.001, clump_p = 1).
Save clumped data to expo_rt_XXX.txt, then process and save to expo_rt_out_XXX.txt.

2.2 Step 2 (reverse MR step2.R)
Set working directory using setwd("").
Read step - 1 - generated exposure data and outcome data for each ID in cathepsin_id.txt.
Harmonize data, conduct MR analysis, calculate statistics, and save results in ID - named directories if MR p - value meets criteria.


3.Multiple - variable MR (multipleMR.R)
Set working directory with setwd("").
Read exposure IDs from cathepsin_id.txt, extract exposure data, read outcome data from the "outcome file", and harmonize data.
Input harmonized data into mr_mvinput. Conduct MR analyses with methods like IVW (fixed and random effects), MR - Egger, and median. Calculate statistics and save results in separate text files.


4.Forest Plot (forestplot.R)
Set working directory using setwd("") where DU.txt is located.
Read data from DU.txt, format p - values and odds ratios with CIs. Use forest function in forestploter to create a forest plot with set aesthetics.
res_fe_ivw_results.txt: IVW (fixed effects) results.
res_re_ivw_results.txt: Another set of IVW (random effects) results.
res_egger_results.txt: MR - Egger results.
res_median_results.txt: Median - based MR results.
3.4 Forest Plot
A visual forest plot representing odds ratios and CIs in DU.txt for result comparison and interpretation.
