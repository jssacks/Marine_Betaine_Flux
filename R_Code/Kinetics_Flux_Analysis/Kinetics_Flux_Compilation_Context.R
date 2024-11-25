##





#load packages
library(tidyverse)





###Define inputs:
NLS.file <- "Intermediates/NLS_Kin_Flux_vals.csv"
WH.file <- "Intermediates/WH_Kin_Flux_vals.csv"
part.conc.file <- "Intermediates/Kinetics_QE_Particulate_Quantification_Results.csv"
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"
match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"
enviro.conc.file <- "Intermediates/Dissolved_Final_Quant_QCed.csv"




# Section 1: Load in datasets --------------------------------------------------------


###Load in Kinetics and Flux data:
nls.dat <- read_csv(NLS.file)
wh.dat <- read_csv(WH.file)



##Load in primary productivity data:
pp.dat <- read_csv(meta.data.file) %>%
  select(KinExp_ID, pp_13c, pp_14c) %>% 
  mutate(pp_13c = replace_na(pp_13c, 0),
         pp_14c = replace_na(pp_14c, 0)) %>%
  mutate(pp.uM.C.perday = case_when(pp_13c > pp_14c ~ pp_13c,
                                    pp_14c > pp_13c ~ pp_14c)) %>%
  select(KinExp_ID, pp.uM.C.perday) %>%
  separate(KinExp_ID, into = c("Cruise", "exp")) %>%
  filter(!exp == "UKG2")



####load in particulate concentration data:
p.conc.dat <- read_csv(part.conc.file) %>%
  mutate(Compound = str_replace(Compound, "Glycine betaine", "GBT")) %>%
  rename("exp" = Exp) %>%
  filter(!exp == "UKG2")



###load in dissolved metabolite data:


#####__Pull in environmental dissolved metabolite data:
diss.dat <- read_csv(enviro.conc.file) %>%
  mutate(Diss_Samp_ID = SampID) 


## organize dissolved data to only include ust samples collected in this study:
exp.enviro.dat <- read_csv(match.file) %>%
  select(-Compound) %>%
  # mutate(Compound = str_replace(Compound, "GBt", "Glycine betaine")) %>%
  left_join(diss.dat) %>%
  # select(KinExp_ID, Diss.Conc.nM) %>%
  group_by(KinExp_ID, Compound) %>%
  reframe(Mean.Diss.Conc.nM = mean(Diss.Conc.nM),
          SD.Diss.Conc.nM = sd(Diss.Conc.nM)) %>%
  separate(KinExp_ID, into=c("Cruise", "exp")) %>%
  filter(!exp == "UKG2")





# Section 2: Contextualize Flux data as % of particulate pool, dissolved pool, and primary production --------------------------------------------------------


#### NLS
nls.context <- nls.dat %>%
  mutate(Compound = case_when(str_detect(exp, "UKG") ~ "GBT",
                              TRUE ~ "Homarine")) %>%
  left_join(., p.conc.dat) %>%
  left_join(., pp.dat) %>%
  mutate(C = case_when(Compound == "Homarine" ~7,
                       Compound == "GBT" ~ 5)) %>%
  mutate(flux_nM_C_day = mm_flux_nM_day*C,
         flux_sd_C = mm_flux_sd_nM_day*C) %>%
  mutate(Flux_Perc_of_Diss = (mm_flux_nM_day/Mean.Diss.Conc.nM)*100,
         Flux_Perc_of_Diss_Error = Flux_Perc_of_Diss*sqrt((mm_flux_sd_nM_day/mm_flux_nM_day)^2+(SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2)) %>%
  mutate(Flux_Perc_of_Part = (mm_flux_nM_day/mean.part.conc.nM)*100,
         Flux_Perc_of_Part_Error = Flux_Perc_of_Part*sqrt((mm_flux_sd_nM_day/mm_flux_nM_day)^2+(sd.part.conc.nM/mean.part.conc.nM)^2)) %>%
  mutate(flux_Perc_of_PP = (flux_nM_C_day/pp.uM.C.perday)*100,
         flux_Perc_of_PP_Error = flux_Perc_of_PP*sqrt((flux_sd_C/flux_nM_C_day)^2))

write_csv(nls.context, file = "Intermediates/NLS_Kin_Flux_context.csv")

#### Wright-Hobbie
wh.context <- wh.dat %>%
  mutate(Compound = case_when(str_detect(exp, "UKG") ~ "GBT",
                              TRUE ~ "Homarine")) %>%
  select(-cruise_exp) %>%
  rename(Cruise = cruise) %>%
  left_join(., p.conc.dat) %>%
  left_join(., pp.dat) %>%
  mutate(C = case_when(Compound == "Homarine" ~7,
                       Compound == "GBT" ~ 5)) %>%
  mutate(flux_nM_C_day = wh_flux_nM_day*C,
         flux_sd_C = wh_flux_sd_nM_day*C) %>%
  mutate(Flux_Perc_of_Diss = (wh_flux_nM_day/Mean.Diss.Conc.nM)*100,
         Flux_Perc_of_Diss_Error = Flux_Perc_of_Diss*sqrt((wh_flux_sd_nM_day/wh_flux_nM_day)^2+(SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2)) %>%
  mutate(Flux_Perc_of_Part = (wh_flux_nM_day/mean.part.conc.nM)*100,
         Flux_Perc_of_Part_Error = Flux_Perc_of_Part*sqrt((wh_flux_sd_nM_day/wh_flux_nM_day)^2+(sd.part.conc.nM/mean.part.conc.nM)^2)) %>%
  mutate(flux_Perc_of_PP = (flux_nM_C_day/pp.uM.C.perday)*100,
         flux_Perc_of_PP_Error = flux_Perc_of_PP*sqrt((flux_sd_C/flux_nM_C_day)^2))

write_csv(wh.context, file = "Intermediates/WH_Kin_Flux_context.csv")




# Section 3: Organize Dissolved data for paper --------------------------------------------------------

###Explore Enviro metab data to summarize for paper:

diss.enviro.sum <- exp.enviro.dat %>%
  group_by(Cruise, exp) %>%
  mutate(Tot.Conc.nM = sum(Mean.Diss.Conc.nM),
         Tot.Conc.nM.SD = sum(SD.Diss.Conc.nM),
         Rel.Abun = Mean.Diss.Conc.nM/Tot.Conc.nM,
         Rel.Abun.SD = Rel.Abun*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(Tot.Conc.nM.SD/Tot.Conc.nM)^2))
  
write_csv(diss.enviro.sum, file = "Intermediates/Dissolved_Metab_Environment_Summary.csv")



# exp.enviro.dat <- read_csv(match.file) %>%
#   select(-Compound) %>%
#   # mutate(Compound = str_replace(Compound, "GBt", "Glycine betaine")) %>%
#   left_join(diss.dat) %>%
#   # select(KinExp_ID, Diss.Conc.nM) %>%
#   group_by(KinExp_ID, Compound) %>%
#   reframe(Mean.Diss.Conc.nM = mean(Diss.Conc.nM),
#           SD.Diss.Conc.nM = sd(Diss.Conc.nM)) %>%
#   separate(KinExp_ID, into=c("Cruise", "exp"))
# 
# write_csv(exp.enviro.dat, file = "Intermediates/Kinetics_Exp_Enviro_Metabolome_Dat.csv")



  
  
###_____ Section 4: Kt/S ratio Calculations: ___________________________________________________
Kt_S_ratio.dat <- nls.context %>%
  select(Cruise, exp, Compound, mean_ks, sd_ks) %>%
  left_join(., diss.enviro.sum %>% 
              select(Cruise, exp, Compound, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, Rel.Abun, Rel.Abun.SD) %>%
              mutate(Compound = str_replace(Compound, "Glycine betaine", "GBT"))) %>% 
  mutate(kt_s_ratio = mean_ks/Mean.Diss.Conc.nM,
         kt_s_ratio_error = kt_s_ratio*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(sd_ks/mean_ks)^2),
         kt_s_ratio_error_fig_max = kt_s_ratio+kt_s_ratio_error,
         kt_s_ratio_error_fig_min = case_when(kt_s_ratio-kt_s_ratio_error < 0 ~ 0,
                                              kt_s_ratio-kt_s_ratio_error >= 0 ~ kt_s_ratio-kt_s_ratio_error)) 

write_csv(Kt_S_ratio.dat, file = "Intermediates/Kt_S_ratio_results.csv")







