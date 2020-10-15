library(gprofiler2)

interesse_regs <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                    "NACC2", "ZNF385B", "NR1D1", "LHX4")

subgroups <- c("g3", "g4", "shh")

load(paste0("../expression_analysis/rdata_files/network/", subgroup, "_rtn.RData"))

regulon_content <- rtna@listOfRegulons[interesse_regs]

for(i in 1:length(regulon_content)){
  write.csv(regulon_content[i], 
            file = paste0("./regulon_content/", 
                          names(regulon_content[i]), ".csv"))
}

terms <- NULL

for(i in 1:length(regulon_content)){
  
  gostres <- gost(query = regulon_content[i],
                  organism = "hsapiens", 
                  sources = c("GO:BP", "KEGG", "REAC"))
  
  
  terms <- append(terms, gostres[["result"]][["term_name"]])
  
  assign(names(regulon_content[i]), gostres)
}


tableau <- sort(table(terms), decreasing = T)[1:20]
