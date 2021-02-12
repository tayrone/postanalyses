library(gprofiler2)
library(gdata)
library(tidyverse)

interesse_regs <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                    "NACC2", "ZNF385B", "NR1D1", "LHX4")

subgroups <- c("g3", "g4", "shh")


#---- Writes down csv files with the content of each regulon, for possible 
# further investigations. Also runs gene ontology analysis for each one of them. 
# Terms are analyzed afterwards. ----

regulon_content <- function(subgroup){

  load(paste0("../expression_analysis/rdata_files/network/", subgroup, 
              "_rtn.RData"))
  
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
                    sources = c("GO:BP", "KEGG", "REAC"),
                    user_threshold = 0.01,
                    correction_method = "bonferroni")
    
    terms <- append(terms, gostres[["result"]][["term_name"]])
    
    assign(names(regulon_content[i]), gostres)
  }
  
  return(unique(terms))

}

terms <- sapply(subgroups, regulon_content)

terms <- as.data.frame(table(stack(terms)))
colnames(terms) <- c("term", "subgroup", "is_present")


#---- As shown below, each term shows up once for each subgroup ----

terms %>% 
  group_by(term) %>% 
  count() %>% 
  filter(n != 4)


#---- Terms that appear in one analysis, only ----

exclusive_terms <- terms %>% 
  group_by(term) %>% 
  mutate(term_freq = sum(is_present)) %>% 
  filter(term_freq == 1 & is_present == 1) %>% 
  select(term, subgroup)


#---- The amount of exclusive terms, by analysis ----

exclusive_terms %>% 
  group_by(subgroup) %>% 
  transmute(n()) %>% 
  distinct()


#---- Terms that appear in all four analyses ----

universal_terms <- terms %>% 
  group_by(term) %>% 
  mutate(term_freq = sum(is_present)) %>% 
  filter(term_freq == 4) %>% 
  select(term) %>% 
  distinct() 

