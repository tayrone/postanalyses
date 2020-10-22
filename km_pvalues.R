library(gdata)

interest_regulons <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                       "NACC2", "ZNF385B", "NR1D1", "LHX4")

subgroups <- c("g3", "g4", "shh", "no_wnt")


km_pvalues <- function(subgroup){
  load(paste0("../expression_analysis/rdata_files/survival/", subgroup, 
              "_survival.RData"))
  
  
  rtns@results$KM$Table[interest_regulons, "Adjusted.Pvalue"]
}

p_values <- lapply(subgroups, km_pvalues)

p_values <- as.data.frame(p_values, col.names = subgroups, row.names = interest_regulons)  
