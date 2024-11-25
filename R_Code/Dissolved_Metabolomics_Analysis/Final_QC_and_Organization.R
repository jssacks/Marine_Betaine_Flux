
library(tidyverse)

#Finalizing Script to organize, apply qc, and combine distribution datasets:

#inputs

#dissolved
d.quant.file <- "Intermediates/Dissolved_Quantified_Data.csv"
d.qc.file.remove <- "Intermediates/Dissolved_osmo_QCdat_samplesremoved.csv"
d.qc.file.impute <- "Intermediates/Dissolved_osmo_QCdat_blanksimputed.csv"
#d.edgecase.file <- "Intermediates/GBT_Quant_for_KM1906_GBT_F2_edge_case.csv"
d.lod.file <- "Intermediates/Dissolved_Blk_LOD_Concentrations.csv"






#Finalize datasets: 



#Dissolved (apply QC and blank subtraction)

#load in datasets
d.q.dat <- read_csv(d.quant.file) %>%
  rename("Rep" = SampID,
         "Compound" = Name)
d.qc.i <- read_csv(d.qc.file.impute)

# d.edgecase.dat <- read_csv(d.edgecase.file) %>%
#   rename("redone.conc" = EE.adjust.conc) %>%
#   select(Rep, Compound, redone.conc) %>%
#   filter(Compound == "Glycine betaine")
# #  mutate(SampID = "KM1906_GBT_F2_T0",
# #         min.area.flag = NA,
# #         blk.lod.flag = NA,
# #         smp.remove = NA


d.lod.dat <- read_csv(d.lod.file) %>%
  rename("Compound" = Name)

##incorporate blk data imputation and blk subtraction into final dataset:
d.q.i.dat <- d.qc.i %>%
  filter(!str_detect(Rep, "Std")) %>%
  filter(!str_detect(Rep, "Blk")) %>%
  filter(!str_detect(Rep, "Poo")) %>%
  left_join(., d.q.dat) %>%
  left_join(., d.lod.dat) %>%
  # left_join(., d.edgecase.dat) %>%             #Add in edgecase dat form weird GBT sample
  # mutate(EE.adjust.conc = case_when(!is.na(redone.conc) ~ redone.conc,
  #                                   TRUE ~ EE.adjust.conc)) %>%
  # select(-redone.conc) %>%
  mutate(Diss.Conc.nM = case_when(blk.lod.flag == "Flag" ~ EE.adjust.Blk.Av/2,
                                  TRUE ~ EE.adjust.conc - EE.adjust.Blk.Av)) %>%
  mutate(Diss.Nmol.C = Diss.Conc.nM*C,
         Diss.Nmol.N = Diss.Conc.nM*N) %>%
  rename("LOD.nM" = EE.adjust.lod) %>%
  mutate(LOD.nM.blk.sub = LOD.nM - EE.adjust.Blk.Av) %>%
  select(Rep, SampID, replicate, Cruise, Compound, min.area.flag, blk.lod.flag, smp.remove,
         Diss.Conc.nM, Diss.Nmol.C, Diss.Nmol.N, LOD.nM, LOD.nM.blk.sub)
         

#write final dataset to csv:
write_csv(d.q.i.dat, file = "Intermediates/Dissolved_Final_Quant_QCed.csv")

