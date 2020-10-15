library(VennDiagram)
library(RTN)
library(gdata)
library(SuperExactTest)


subgroups <- c("g3", "g4", "shh")
names(subgroups) <- c("Group3", "Group4", "SHH")

for(i in 1:length(subgroups)){
  load(paste0("../expression_analysis/rdata_files/network/", subgroups[i], 
              "_rtn.RData"))
  
  mrs <- tna.get(rtna, what = "mra")
  
  load(paste0("../methylation/rdata_files/", names(subgroups)[i], 
              "/4_dm_regulons.RData"))
  
  gdata::keep(mrs, hm_regulons, subgroups, i, sure = T)
  
  load(paste0("../expression_analysis/rdata_files/survival/", subgroups[i], 
              "_survival.RData"))
  
  gdata::keep(mrs, hm_regulons, hazardous_regulons, subgroups, i, sure = T)
  
  venn.diagram(
    x = list(mrs$Regulon, rownames(hm_regulons), hazardous_regulons),
    category.names = c("MRs" , "DM Regulons" , "Hazardous"),
    main = names(subgroups[i]),
    filename = paste0(subgroups[i], "_venn_diagramm.png"),
    output = TRUE, height = 5000, width = 6500, resolution = 750
  )

}


gdata::keep(subgroups, sure = T)

#----

for(subgroup in subgroups){
  load(paste0("../expression_analysis/rdata_files/network/", subgroup, 
              "_rtn.RData"))
  
 assign(paste0(subgroup, "_mrs"), tna.get(rtna, what = "mra"))
  
}

venn.diagram(
  x = list(g3_mrs$Regulon, g4_mrs$Regulon, shh_mrs$Regulon),
  category.names = c("Grupo 3" , "Grupo 4" , "SHH"),
  main = "Reguladores Mestres",
  filename = "rms_venn_diagramm.png",
  output = TRUE, height = 5000, width = 6500, resolution = 750
)

cat("\n\nMRs: in g3 g4, out shh\n", noquote(
  setdiff(intersect(g3_mrs$Regulon, g4_mrs$Regulon), shh_mrs$Regulon)
  ))

cat("\n\nMRs (g3, g4, shh)\n", noquote(
  intersect(intersect(g3_mrs$Regulon, g4_mrs$Regulon), shh_mrs$Regulon)
))

gdata::keep(subgroups, g3_mrs, g4_mrs, shh_mrs, sure = T)

#----


for(subgroup in names(subgroups)){
  load(paste0("../methylation/rdata_files/", subgroup, 
              "/4_dm_regulons.RData"))
  
  assign(paste0(subgroup, "_dm_regulons"), rownames(hm_regulons))
  
}

venn.diagram(
  x = list(Group3_dm_regulons, Group4_dm_regulons, SHH_dm_regulons),
  category.names = c("Grupo 3" , "Grupo 4" , "SHH"),
  main = "Regulons Diferencialmente Metilados",
  filename = "dm_regulons_venn_diagramm.png",
  output = TRUE, height = 5000, width = 6500, resolution = 750
)

gdata::keep(subgroups, g3_mrs, g4_mrs, shh_mrs, 
            Group3_dm_regulons, Group4_dm_regulons, SHH_dm_regulons,
            sure = T)


#----

g3_interest <- intersect(g3_mrs$Regulon, Group3_dm_regulons)
g4_interest <- intersect(g4_mrs$Regulon, Group4_dm_regulons)
shh_interest <- intersect(shh_mrs$Regulon, SHH_dm_regulons)

cat("\n\nInterest: g3, g4, shh\n",
    noquote(intersect(intersect(g3_interest, g4_interest), shh_interest)))


venn.diagram(
  x = list(g3_interest, g4_interest, shh_interest),
  category.names = c("Grupo 3" , 
                     "Grupo 4" , 
                     "SHH"),
  main = "Regulons de Interesse",
  sub = "Regulons diferencialmente expressos e metilados, simultaneamente",
  filename = "regs_de_interesse_venn_diagramm.png",
  output = TRUE, height = 5000, width = 6500, resolution = 750
)


#gdata::keep(subgroups, sure = T)

load("../expression_analysis/rdata_files/tfs.RData")

supertest(list(g3_interest, g4_interest, shh_interest), n = length(tfs))$P.value
