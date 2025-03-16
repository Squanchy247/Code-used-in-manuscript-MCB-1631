install.packages("forestploter")
install.packages("grid")

library(grid)
library(forestploter)

setwd("")

result = read.table("DU.txt", header = T, sep = "\t")
result$pval = ifelse(result$pval < 0.001, "<0.001", sprintf("%.4f", result$pval))
result$`OR (95% CI)` = ifelse(is.na(result$or), "", sprintf("%.4f (%.4f - %.4f)",
                                                            result$or, result$or_lci95, 
                                                            result$or_uci95))
result$'P-value' = result$pval

forest(result[, c(1, 6:8)],
       est = as.numeric(result$or),
       lower = as.numeric(result$or_lci95),
       upper = as.numeric(result$or_uci95),
       sizes = 0.3,
       ci_column = 2,
       ref_line = 1,
       xlim = c(0.05, 2),
       theme = forest_theme(ci_Theight = 0.2))