



library(tidyverse)
library(scales)




#define inputs
nls.file <- "Intermediates/NLS_Kin_Flux_context.csv"
wh.file <- "Intermediates/WH_Kin_Flux_context.csv"
#f.rel.file <- "Intermediates/flux_contextualizaing_data"
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"
Kt_S_ratio.file <- "Intermediates/Kt_S_ratio_results.csv"
enviro.metab.file <- "Intermediates/Dissolved_Metab_Environment_Summary.csv"
#dnorm.file <-  "Intermediates/DNorm_dat.csv"



###Make Table 1 (Kinetics Parameters)

#Get kinetics data
nls.dat <- read_csv(nls.file)

#Get Region dat:
region.dat <- read_csv(meta.data.file) %>%
  select(KinExp_ID, Region) %>%
  separate(KinExp_ID, into = c("Cruise", "exp"), remove = FALSE) %>%
  mutate(Region = case_when(Region == "Gyre" ~ "NPSG",
                            Region == "Coastal" ~ "PS",
                            Region == "Equatorial" ~ "Eq",
                            TRUE ~ Region)) %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("NPTZ", "NPSG", "Eq", "PS"))) 


#Compile data:
kin.table <- nls.dat %>%
  left_join(., region.dat) %>%
  select(Region, Cruise, exp, Compound, mean.part.conc.nM, sd.part.conc.nM, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, mean_ks, sd_ks, mean_vmax, sd_vamx) %>%
  rename("sd_vmax" = sd_vamx) %>%
  mutate(Mean.Part.Conc.nM = print(formatC(signif(round(mean.part.conc.nM, 3), digits=3), digits=3,format="fg")),
         SD.Part.Conc.nM = print(formatC(signif(round(sd.part.conc.nM,3), digits=2), digits=2,format="fg")),
         Mean.Diss.Conc.nM = print(formatC(signif(round(Mean.Diss.Conc.nM,3), digits=3), digits=3,format="fg")),
         SD.Diss.Conc.nM = print(formatC(signif(round(SD.Diss.Conc.nM,3), digits=2), digits=2,format="fg")),
         mean_ks = print(formatC(signif(mean_ks, digits=3), digits=3,format="fg")),
         sd_ks = print(formatC(signif(sd_ks, digits=2), digits=2,format="fg")),
         mean_vmax = print(formatC(signif(round(mean_vmax,3), digits=3), digits=3,format="fg")),
         sd_vmax = print(formatC(signif(round(sd_vmax,3), digits=2), digits=2,format="fg"))) %>%
  mutate(P.Conc.nM = paste(Mean.Part.Conc.nM, SD.Part.Conc.nM, sep = "\u00B1"),
         D.Conc.nM = paste(Mean.Diss.Conc.nM, SD.Diss.Conc.nM, sep = "\u00B1"),
         ks = paste(mean_ks, sd_ks, sep = "\u00B1"),
         vmax = paste(mean_vmax, sd_vmax, sep = "\u00B1")) %>%
  select(Region, Cruise, exp, Compound, P.Conc.nM, D.Conc.nM, ks, vmax) %>%
  rename("Experiment" = exp,
         "Dissolved Concentration (nM)" = D.Conc.nM,
         "Particulate Concentration (nM)" = P.Conc.nM) %>%
  rename("Kt (nM)" = ks,
         "Vmax nmol L^-1 day^-1" = vmax) %>%
  arrange(Compound)

#export
write_csv(kin.table, file = "Tables/Main_Text_Table1_Kinetics.csv")




###TT and FLux Table (Main Text Table 2)
flux.table <- nls.dat %>%
  left_join(., region.dat) %>%
  select(Region, Cruise, exp, Compound, mm_tt, mm_tt_sd, mm_flux_nM_day, mm_flux_sd_nM_day) %>%
  mutate(carbon = case_when(Compound == "GBT" ~ 5,
                            Compound == "Homarine" ~ 7),
         C_flux_nM_day = mm_flux_nM_day*carbon,
         C_flux_nM_day_sd = mm_flux_sd_nM_day*carbon) %>%
  mutate(TT = print(formatC(signif(mm_tt, digits=3), digits=3,format="fg")),
         TT_sd = print(formatC(signif(mm_tt_sd, digits=2), digits=2,format="fg")),
         flux = print(formatC(signif(round(mm_flux_nM_day,3), digits=3), digits=3,format="fg")),
         flux_sd = print(formatC(signif(round(mm_flux_sd_nM_day,3), digits=2), digits=2,format="fg")),
         c_flux = print(formatC(signif(round(C_flux_nM_day,3), digits=3), digits=3,format="fg")),
         c_flux_sd = print(formatC(signif(round(C_flux_nM_day_sd,3), digits=2), digits=2,format="fg"))) %>%
  mutate(TT = paste(TT, TT_sd, sep = "\u00B1"),
         Flux = paste(flux, flux_sd, sep = "\u00B1"),
         C_Flux = paste(c_flux, c_flux_sd, sep = "\u00B1")) %>%
  select(Region, Cruise, exp, Compound, TT, Flux, C_Flux) %>%
  rename("Experiment" = exp,
         "Turnover Time (h)" = TT,
         "Flux (nmol compound L^-1 day^-1)" = Flux,
         "Carbon Flux (nmol C L^-1 day^-1)" = C_Flux) %>%
  arrange(Compound)

#export
write_csv(flux.table, file = "Tables/Main_Text_Table2_Fluxes.csv")






###Make Wright-Hobbie TT Table (Supplemental):
WH.table <- read_csv(wh.file) %>%
  select(Cruise, exp, Compound, wh_mean_kt_sn, wh_sd_kt_sn, wh_mean_vmax, wh_sd_vmax, wh_mean_tt, wh_sd_tt) %>%
  mutate(mean_ktsn = print(formatC(signif(wh_mean_kt_sn, digits=3), digits=3,format="fg")),
         sd_ktsn = print(formatC(signif(wh_sd_kt_sn, digits=2), digits=2,format="fg")),
         mean_vmax = print(formatC(signif(wh_mean_vmax, digits=3), digits=3,format="fg")),
         sd_vmax = print(formatC(signif(wh_sd_vmax, digits=2), digits=2,format="fg")),
         TT = print(formatC(signif(wh_mean_tt, digits=3), digits=3,format="fg")),
         TT_sd = print(formatC(signif(wh_sd_tt, digits=2), digits=2,format="fg"))) %>%
  select(Cruise, exp, Compound, mean_ktsn, sd_ktsn, mean_vmax, sd_vmax, TT, TT_sd) %>%
  rename("Experiment" = exp,
         "Kt+Sn (nM)" = mean_ktsn,
         "Kt+Sn SD (nM)" = sd_ktsn,
         "Vmax (nmol L^-1 day^-1)" = mean_vmax,
         "Vmax SD (nmol L^-1 day^-1)" = sd_vmax,
         "TT (h)" = TT,
         "TT SD (h)" = TT_sd) %>%
  arrange(Compound) 

#export
write_csv(WH.table, file = "Tables/WH_Supplemental_Table.csv")






###Make Flux in Context Table (Percent of Primary Production, etc.): 
f.rel.table <- read_csv(nls.file) %>%
  select(Cruise, exp, Compound, mean.part.conc.nM, sd.part.conc.nM, Flux_Perc_of_Diss, Flux_Perc_of_Diss_Error,
         Flux_Perc_of_Part, Flux_Perc_of_Part_Error, pp.uM.C.perday, flux_Perc_of_PP, flux_Perc_of_PP_Error) %>%
  mutate(mean.part.conc.nM = print(formatC(signif(mean.part.conc.nM, digits=3), digits=3,format="fg")),
         sd.part.conc.nM = print(formatC(signif(sd.part.conc.nM, digits=2), digits=2,format="fg")),
         Flux_Perc_of_Diss = print(formatC(signif(Flux_Perc_of_Diss, digits=3), digits=3,format="fg")),
         Flux_Perc_of_Diss_Error = print(formatC(signif(Flux_Perc_of_Diss_Error, digits=2), digits=2,format="fg")),
         Flux_Perc_of_Part = print(formatC(signif(Flux_Perc_of_Part, digits=3), digits=3,format="fg")),
         Flux_Perc_of_Part_Error = print(formatC(signif(Flux_Perc_of_Part_Error, digits=2), digits=2,format="fg")),
         pp.uM.C.perday = print(formatC(signif(pp.uM.C.perday, digits=3), digits=3,format="fg")),
         flux_Perc_of_PP = print(formatC(signif(flux_Perc_of_PP, digits=3), digits=3,format="fg")),
         flux_Perc_of_PP_Error = print(formatC(signif(flux_Perc_of_PP_Error, digits=2), digits=2,format="fg",))) %>%
  rename("Experiment" = exp,
         "Average Particulate Concentration (nM)" = mean.part.conc.nM,
         "Standard Deviation of Particulate Concentration (nM)" = sd.part.conc.nM, 
         "Daily Flux as Percentage of Dissolved Pool (%)" = Flux_Perc_of_Diss,
         "Daily Flux as Percentage of Dissolved Pool Error (%)" = Flux_Perc_of_Diss_Error,
         "Daily Flux as Percentage of Particulate Pool (%)" = Flux_Perc_of_Part,
         "Daily Flux as Percentage of Particulate Pool Error (%)" = Flux_Perc_of_Part_Error,
         "Estimated Primary Production (uM C day^-1)" = pp.uM.C.perday,
         "Daily Flux as Percentage of Primary Production (%)" = flux_Perc_of_PP,
         "Daily Flux as Percentage of Primary Production Error (%)" = flux_Perc_of_PP_Error) %>%
  arrange(Compound)

#export:
write_csv(f.rel.table, file = "Tables/Flux_in_Context_Supplemental_Table.csv")



###Make Kt_S_ratio file:

#organize and compile dissolved data and kt_s_ratio data
kt_s_ratio.dat <- read_csv(Kt_S_ratio.file) %>%
  select(Cruise, exp, Compound, kt_s_ratio, kt_s_ratio_error)

dissolved.metab.dat <- read_csv(enviro.metab.file) %>%
  select(Cruise, exp, Compound, Mean.Diss.Conc.nM, SD.Diss.Conc.nM) %>%
  filter(!Compound == "Trimethylamine N-oxide") %>%
  mutate(Compound = str_replace(Compound, "Glycine betaine", "GBT")) %>% 
  group_by(Cruise, exp) %>%
  mutate(Tot.Sulf.Bet.Conc.nM = sum(Mean.Diss.Conc.nM),
         Tot.Sulf.Bet.Conc.nM.SD = sum(SD.Diss.Conc.nM),
         Perc.tot.pool = Mean.Diss.Conc.nM/Tot.Sulf.Bet.Conc.nM*100,
         Perc.tot.pool.SD =  Perc.tot.pool*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(Tot.Sulf.Bet.Conc.nM.SD/Tot.Sulf.Bet.Conc.nM)^2))

kt_s_ratio.diss.dat <- left_join(kt_s_ratio.dat, dissolved.metab.dat)

####Make Kt S ratio table:
kt_s_ratio.table <- kt_s_ratio.diss.dat %>%
  select(Cruise, exp, Compound, kt_s_ratio, kt_s_ratio_error, Tot.Sulf.Bet.Conc.nM, Tot.Sulf.Bet.Conc.nM.SD, Perc.tot.pool, Perc.tot.pool.SD) %>%
  mutate(kt_s_ratio = print(formatC(signif(kt_s_ratio, digits=3), digits=3,format="fg")),
         kt_s_ratio_error = print(formatC(signif(kt_s_ratio_error, digits=2), digits=3,format="fg")),
         Tot.Sulf.Bet.Conc.nM = print(formatC(signif(Tot.Sulf.Bet.Conc.nM, digits=3), digits=3,format="fg")),
         Tot.Sulf.Bet.Conc.nM.SD = print(formatC(signif(Tot.Sulf.Bet.Conc.nM.SD, digits=2), digits=2,format="fg")),
         Perc.tot.pool = print(formatC(signif(Perc.tot.pool, digits=3), digits=3,format="fg")),
         Perc.tot.pool.SD = print(formatC(signif(Perc.tot.pool.SD, digits=2), digits=2,format="fg"))) %>%
  rename("Experiment" = exp,
         "Kt/Sn Ratio" = kt_s_ratio,
         "Kt/Sn Ratio Error" = kt_s_ratio_error,
         "Total Concentration of Sulfoniums and Betaines (nM)" = Tot.Sulf.Bet.Conc.nM,
         "Total Concentration of Sulfoniums and Betaines Error (nM)" = Tot.Sulf.Bet.Conc.nM.SD,
         "Compound Percentage of Total Pool (%)" = Perc.tot.pool,
         "Compound Percentage of Total Pool Error (%)" = Perc.tot.pool.SD) %>%
  arrange(Compound)

write_csv(kt_s_ratio.table, file = "Tables/Kt_Sn_Ratio_Table.csv")








