library(gdata)
library(tidyverse)
library(RTNsurvival)

interest_regulons <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                       "NACC2", "ZNF385B", "NR1D1", "LHX4")

subgroups <- c("g3", "g4", "shh", "no_wnt")


#---- Defines Kaplan-Meier estimation p-value (aka log-rank test p-value), 
# for all regulons of interest ----

km_pvalues <- function(subgroup){
  
  load(paste0("../expression_analysis/rdata_files/survival/", subgroup, 
              "_survival.RData"))
  
  rtns@results$KM$Table[interest_regulons, "Adjusted.Pvalue"]
}

p_values <- lapply(subgroups, km_pvalues)

p_values <- as.data.frame(p_values, col.names = subgroups, row.names = interest_regulons)  


#---- Plots an histogram for log rank test adjusted p-values ----

subgroups <- c("g3", "g4", "shh", "no_wnt")

km_pvalues <- function(subgroup){
  load(paste0("../expression_analysis/rdata_files/survival/", subgroup, 
              "_survival.RData"))
  
  rtns@results[["KM"]][["Table"]] %>% 
    select(Regulons, Adjusted.Pvalue) %>% 
    ggplot() +
    geom_histogram(aes(Adjusted.Pvalue), binwidth = 0.1, center = 0.05) +
    theme_minimal() +
    labs(title = paste0("Subgroup: ", subgroup), y = "Number of regulons")
  
    ggsave(file = paste0("./histogram_plots/", subgroup,".png"))
}

pvalue_histograms <- lapply(subgroups, km_pvalues)

#gdata::keep(pvalue_histograms, interest_regulons, sure = T)

