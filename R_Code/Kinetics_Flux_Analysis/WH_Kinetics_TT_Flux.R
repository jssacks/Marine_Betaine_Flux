





#packages
library(tidyverse)
library(lubridate)
library(broom)

source("R_Code/Kinetics_Flux_Analysis/Functions.R")

###inputs
rate.file <- "Intermediates/quantified_uptake_rates.csv"
boysen.rate.file <- "Meta_Data/Data_From_Other_Studies/Boysen2022_GBT_uptake_rates.csv"
match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"
enviro.conc.file <- "Intermediates/Dissolved_Final_Quant_QCed.csv"
#kin.file <- "Intermediates/Kin_Analysis_Output.csv"
#exp.diss.match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"



# Pull in and Organize Data -----------------------------------------------

###Pull in dissolved data
diss.dat <- read_csv(enviro.conc.file) %>%
  mutate(Diss_Samp_ID = SampID) 

exp.enviro.dat <- read_csv(match.file) %>%
  mutate(Compound = str_replace(Compound, "GBt", "Glycine betaine")) %>%
  left_join(diss.dat) %>%
  filter(!Rep == "220602_Smp_KM1906_GBT_F2_T0_A") %>%
  select(KinExp_ID, Diss.Conc.nM) %>%
  group_by(KinExp_ID) %>%
  reframe(Mean.Diss.Conc.nM = mean(Diss.Conc.nM),
          SD.Diss.Conc.nM = sd(Diss.Conc.nM)) %>%
  separate(KinExp_ID, into=c("cruise", "exp"), remove = FALSE) %>%
  rename("cruise_exp" = KinExp_ID) %>%
  filter(!exp == "UKG2")


#pull in uptake experiment rate data:
dat.rate = read_csv(rate.file) 

###Pull in GBT uptake rate data from Gradients 3 cruise from Boysen et al. 2022
boysen.dat <- read_csv(boysen.rate.file) %>%
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
  select(-Replicate.Name) %>%
  filter(!exp == "UKG2")


##Subset data and combine boysen et al data with data from this study
dat.uptake <- dat.rate %>%
  filter(!exp == "UCH1") %>%
  filter(!str_detect(SampID, "Blk")) %>%
  filter(!str_detect(SampID, "blk")) %>%
  select(cruise, exp, treatment_conc, rep, Fragment_mz, nM.per.hour.1, se.nM.per.hour.1) %>%
  rbind(boysen.dat)




# Apply Wright-Hobbie Transformation to Calculate Kinetics Parameters (with Monte Carlo Permutation) --------

#define number of models
reps <- 1000
set.seed(981)

#generate MC dataset for uptake kinetics experiments, remove RC104 samples >1000 nM because of high variability
dat.mc <- dat.uptake %>%
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
  mutate(time_over_fract = 1/(mod.val/treatment_conc)) %>%
  mutate(remove = case_when(cruise == "RC104" & treatment_conc > 1000 ~ "Yes",
                            TRUE ~ "No")) %>%
  filter(remove == "No") %>%
  select(-remove) %>%
  group_by(cruise, exp, model)


#Apply WH transformation and retrieve slope and intercept
wh.calc <- do(dat.mc,
              tidy(
                lm(time_over_fract ~ treatment_conc, data = .)))


##Remove outlier models 
wh.qc <- wh.calc %>%
  unite(c(exp, cruise), col = "Exp_Name", remove = FALSE) %>%
  group_by(exp, cruise) %>%
  mutate(mean.std.error = mean(std.error),
         sd.std.error = sd(std.error)) %>%
  group_by(exp, cruise, model) %>%
  mutate(se.flag.high = std.error > (mean.std.error+(3*sd.std.error)),
         se.flag.low = std.error < (mean.std.error - (3*sd.std.error))) %>%
  filter(!any(se.flag.high == TRUE),
         !any(se.flag.low == TRUE)) %>%
  ungroup()


#Calculate Kinetics Values and Turnover Times 
wh.kin.vals <- wh.qc %>%
  select(Exp_Name, exp, cruise, model, term, estimate) %>%
  pivot_wider(id_cols = c("Exp_Name", "exp", "cruise", "model"), names_from = term, values_from = estimate) %>%
  rename("intercept" = `(Intercept)`,
         "slope" = treatment_conc) %>%
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


#### Select just final kinetics parameters:
wh.kin.sml <- wh.kin.vals %>%
  select(cruise, exp, wh_mean_tt, wh_sd_tt, wh_mean_vmax, wh_sd_vmax, wh_mean_kt_sn, wh_sd_kt_sn)


# Combine Wright-Hobbie Kinetics with Dissolved Concentrations to Calculate Fluxes --------

###Wright-Hobbie Flux calculations:
wh.flux.dat <- left_join(wh.kin.sml, exp.enviro.dat) %>%
  mutate(Compound = case_when(str_detect(exp, "UKH") ~ "Homarine",
                              TRUE ~ "GBT")) %>%
  mutate(wh_flux_nM_day = Mean.Diss.Conc.nM/(wh_mean_tt/24),
         wh_flux_sd_nM_day = wh_flux_nM_day*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(wh_sd_tt/wh_mean_tt)^2),
         wh_rel_flux_error = wh_flux_sd_nM_day/wh_flux_nM_day) %>%
  unite(c(cruise, exp), col = cruise_exp, remove = FALSE) %>%
  mutate(Compound = case_when(str_detect(cruise_exp, "UKH") ~ "Homarine",
                              TRUE ~ "GBT")) %>%
  select(cruise_exp, cruise, exp, Compound, everything())


#Export WH Kinetics and Flux Values
write_csv(wh.flux.dat, "Intermediates/WH_Kin_Flux_vals.csv")
