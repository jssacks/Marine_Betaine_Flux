
#packages
library(tidyverse)
library(lubridate)
library(broom)

source("R_Code/Kinetics_DataProcessing/Functions.R")

###inputs
rate.file <- "Intermediates/quantified_uptake_rates.csv"
match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"
enviro.conc.file <- "Intermediates/Dissolved_Final_Quant_QCed.csv"
kin.file <- "Intermediates/Kin_Analysis_Output.csv"
#exp.diss.match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"

###Pull in dissolved data
exp.enviro.dat <- read_csv(match.file) %>%
  mutate(Compound = str_replace(Compound, "GBt", "Glycine betaine")) %>%
  left_join(diss.dat) %>%
  filter(!Rep == "220602_Smp_KM1906_GBT_F2_T0_A") %>%
  select(KinExp_ID, Diss.Conc.nM) %>%
  group_by(KinExp_ID) %>%
  reframe(Mean.Diss.Conc.nM = mean(Diss.Conc.nM),
          SD.Diss.Conc.nM = sd(Diss.Conc.nM)) %>%
  separate(KinExp_ID, into=c("cruise", "exp"), remove = FALSE) %>%
  rename("cruise_exp" = KinExp_ID)







#pull in uptake experiment rate data:
dat.rate = read_csv(rate.file) 

###Pull in GBT uptake rate data from Gradients 3 cruise from Boysen et al. 2022
boysen.dat <- read_csv("Meta_Data/Data_From_Other_Studies/Boysen2022_GBT_uptake_rates.csv") %>%
  select(Replicate.Name, Treatment, replicate, correct.Prediction.nmoles.per.L.per.hr, SE.correct.rate) %>%
  rename("treatment_conc" = Treatment,
         "nM.per.hour.1" = correct.Prediction.nmoles.per.L.per.hr,
         "se.nM.per.hour.1" = SE.correct.rate) %>%
  mutate(cruise = "KM1906") %>%
  mutate(exp = case_when(str_detect(Replicate.Name, "K1") ~ "UKG1",
                         str_detect(Replicate.Name, "K2") ~ "UKG2")) %>%
  filter(replicate %in% c("A", "B", "C")) %>%
  rename("rep" = replicate) %>%
  mutate("Fragment_mz" = 62.8977) %>%
  select(-Replicate.Name)



#_______Calculate Michelis Mentin Uptake Kinetic Parameters by Fitting Nonlinear Equations to Uptake Curves

##Subset data and combine boysen et al data with data from this study
dat.MM <- dat.rate %>%
  filter(!exp == "UCH1") %>%
  filter(!str_detect(SampID, "Blk")) %>%
  filter(!str_detect(SampID, "blk")) %>%
  select(cruise, exp, treatment_conc, rep, Fragment_mz, nM.per.hour.1, se.nM.per.hour.1) %>%
  rbind(boysen.dat)

#define number of models
reps <- 1000
set.seed(981)

#generate MC dataset for uptake kinetics experiments, remove RC104 samples >1000 nM because of high variability
dat.tt.mc <- dat.MM %>%
  filter(Fragment_mz %in% c(97.0954, 62.8977)) %>%
  filter(!treatment_conc == 0) %>%
  group_by(cruise, exp, nM.per.hour.1, se.nM.per.hour.1) %>%
  mutate(mod.val = list(rnorm(reps, mean = nM.per.hour.1,
                              sd = se.nM.per.hour.1 * sqrt(7)))) %>%
  unnest(mod.val) %>%
  mutate(mod.val = case_when(mod.val < 0 ~ nM.per.hour.1,
                             TRUE ~ mod.val)) %>%
  mutate(model = paste("model",rep(1:reps),sep="_")) %>%
  ungroup() %>%
  mutate(time_over_fract = 1/(mod.val/treatment_conc),
         log10.tof = log10(time_over_fract),
         log10.t.nM = log10(treatment_conc)) %>%
  mutate(log10.tof = time_over_fract,
        log10.t.nM = treatment_conc) %>%
  mutate(remove = case_when(cruise == "RC104" & treatment_conc > 1000 ~ "Yes",
                            TRUE ~ "No")) %>%
  filter(remove == "No") %>%
  select(-remove) %>%
  group_by(cruise, exp, model)



tt.calc <- do(dat.tt.mc,
              tidy(
                lm(log10.tof ~ log10.t.nM, data = .)))

wh.qc <- tt.calc %>%
#  filter(term == "(Intercept)") %>%
#  mutate(tt = 10^estimate) %>%
  unite(c(exp, cruise), col = "Exp_Name", remove = FALSE) %>%
  group_by(exp, cruise) %>%
  mutate(mean.std.error = mean(std.error),
         sd.std.error = sd(std.error)) %>%
  group_by(exp, cruise, model) %>%
  mutate(se.flag.high = std.error > (mean.std.error+(3*sd.std.error)),
         se.flag.low = std.error < (mean.std.error - (3*sd.std.error))) %>%
#  ungroup()# %>%
  filter(!any(se.flag.high == TRUE),
         !any(se.flag.low == TRUE))


wh.kin.vals <- wh.qc %>%
  select(Exp_Name, exp, cruise, model, term, estimate) %>%
#  mutate(estimate = 10^estimate) %>%
  pivot_wider(id_cols = c("Exp_Name", "exp", "cruise", "model"), names_from = term, values_from = estimate) %>%
  rename("intercept" = `(Intercept)`,
         "slope" = log10.t.nM) %>%
  mutate(tt = intercept,
         vmax = 1/slope,
         kt_sn = abs(-intercept/slope)) %>%
  group_by(cruise, exp) %>%
  reframe(wh_mean_tt = mean(tt),
          wh_median_tt = median(tt),
          wh_sd_tt = sd(tt),
          wh_max_tt = max(tt),
          wh_min_tt = min(tt),
          wh_mean_vmax = mean(vmax),
          wh_median_vmax = median(vmax),
          wh_sd_vmax = sd(vmax),
          wh_max_vmax = max(vmax),
          wh_min_vmax = min(vmax),
          wh_mean_kt_sn = mean(kt_sn),
          wh_median_kt_sn = median(kt_sn),
          wh_sd_kt_sn = sd(kt_sn),
          wh_max_kt_sn = max(kt_sn),
          wh_min_kt_sn = min(kt_sn))




wh.kin.final <- wh.kin.vals %>%
  select(cruise, exp, wh_mean_tt, wh_sd_tt, wh_mean_vmax, wh_sd_vmax, wh_mean_kt_sn, wh_sd_kt_sn)

###Wright-Hobbie Flux calculations:
wh.flux.dat <- left_join(wh.kin.final, exp.enviro.dat) %>%
  mutate(Compound = case_when(str_detect(exp, "UKH") ~ "Homarine",
                              TRUE ~ "GBT")) %>%
  mutate(wh_flux_nM_day = Mean.Diss.Conc.nM/(wh_mean_tt/24),
         wh_flux_sd_nM_day = wh_flux_nM_day*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(wh_sd_tt/wh_mean_tt)^2),
         wh_rel_flux_error = wh_flux_sd_nM_day/wh_flux_nM_day) %>%
  unite(c(cruise, exp), col = cruise_exp, remove = FALSE) %>%
  mutate(Compound = case_when(str_detect(cruise_exp, "UKH") ~ "Homarine",
                              TRUE ~ "GBT"))





######____Use Kinetics results and environmental concentrations to estimate TT from MM curves and dissolved conc.

kin.dat <- read_csv(kin.file) %>%
  rename("cruise" = Cruise)

#####Pull in environmental data:
diss.dat <- read_csv(enviro.conc.file) %>%
  mutate(Diss_Samp_ID = SampID) 



tt.mm.calc.dat <- left_join(exp.enviro.dat, kin.dat) %>%
  select(cruise, exp, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, mean_ks, sd_ks, mean_vmax, sd_vamx) %>%
  mutate(v_insitu = (Mean.Diss.Conc.nM*mean_vmax)/(Mean.Diss.Conc.nM+mean_ks),
         v_insitu_sd = v_insitu*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(sd_vamx/mean_vmax)^2+(sd_ks/mean_ks)^2),
         mm_tt = Mean.Diss.Conc.nM/v_insitu,
         mm_tt_sd = mm_tt*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(v_insitu_sd/v_insitu)^2),
         mm_flux = v_insitu,
         mm_flux_sd = v_insitu_sd) %>%
         #  mm_flux = Mean.Diss.Conc.nM/mm_tt,
       #  mm_flux_sd = mm_flux*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(mm_tt_sd/mm_tt)^2)) %>%
  mutate(mm_flux_nM_day = mm_flux*24,
         mm_flux_sd_nM_day = mm_flux_sd*24)



######Combine fluxes output:
all.flux.output <- left_join(tt.mm.calc.dat, wh.flux.dat) %>%
  select(cruise_exp, cruise, exp, Compound, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, v_insitu, v_insitu_sd, 
         mm_tt, mm_tt_sd, mm_flux_nM_day, mm_flux_sd_nM_day, wh_mean_tt, wh_sd_tt, wh_flux_nM_day, wh_flux_sd_nM_day)

###Write Flux data to csv
write_csv(all.flux.output, file = "Intermediates/Flux_Analysis_Output.csv")






##Flux viz
#ggplot(flux.dat, aes(x = cruise_exp, y = flux_nM_day)) +
#  geom_point(size = 3) +
#  geom_errorbar(aes(ymin = flux_nM_day-flux_error, ymax = flux_nM_day+flux_error), width = 0.5) +
#  facet_grid(.~Compound, scales = "free_x", space = "free") +
#  scale_y_log10() +  
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 70, vjust = 0.5)) 
 # scale_y_log10() +


#ggplot(flux.dat, aes(x = cruise_exp, y = mean_tt)) +
#  geom_col(width = 0.6, fill = "gray") +
#  facet_grid(.~Compound, scales = "free_x", space = "free") +
#  theme_bw() +
#  geom_errorbar(aes(ymin = mean_tt-sd_tt, ymax = mean_tt+sd_tt), width = 0.3) +
#  theme(axis.text.x = element_text(angle = 70, vjust = 0.5)) 
  
#  geom_point(size = 3) +
#  geom_errorbar(aes(ymin = flux_nM_day-flux_error, ymax = flux_nM_day+flux_error), width = 0.5) +
#  facet_grid(.~Compound, scales = "free_x", space = "free") +
#  scale_y_log10() +  
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 70, vjust = 0.5)) 








