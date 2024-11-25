



#packages
library(tidyverse)
library(lubridate)
library(broom)

source("R_Code/Kinetics_Flux_Analysis/Functions.R")

###inputs
rate.file <- "Intermediates/quantified_uptake_rates.csv"
match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"
enviro.conc.file <- "Intermediates/Dissolved_Final_Quant_QCed.csv"
#exp.diss.match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"


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
  select(-Replicate.Name) %>%
  filter(!exp == "UKG2")







#_______Calculate Michelis Mentin Uptake Kinetic Parameters by Fitting Nonlinear Equations to Uptake Curves

##Subset data and combine boysen et al data with data from this study
dat.MM <- dat.rate %>%
  filter(!exp == "UCH1") %>%
  filter(!str_detect(SampID, "Blk")) %>%
  filter(!str_detect(SampID, "blk")) %>%
  select(cruise, exp, treatment_conc, rep, Fragment_mz, nM.per.hour.1, se.nM.per.hour.1) %>%
  rbind(boysen.dat)

write_csv(dat.MM, "Intermediates/All_Kin_Exp_Dat.csv")


#_____Visualize and estimate starting parameters for MM models_________________ 

###Visualize data to identify starting parameters 
dat.MM.prelimviz <- dat.MM %>%
  group_by(cruise, exp) %>%
  mutate(max.rate = max(nM.per.hour.1, na.rm = TRUE))


###Examine raw data to define starting predictions for Vmax (a) and Ks (b) 
# the red line is the max uptake rate and likely represents a good prediction for Vmax,
ggplot(dat.MM.prelimviz, aes(x = treatment_conc, y = nM.per.hour.1)) +
  geom_point() +
  facet_wrap(cruise~exp, scales = "free") +
  geom_hline(aes(yintercept = max.rate), color = "red") 

##Inspect plots and manually enter predicted values into dataframe, 
# enter values in order of cruise, (cruise), experiment (exp) for a_pred (Vmax prediction) and b_pred (Ks prediction)
model.starting.estiamtes <- data.frame(
  cruise = c("KM1906", "RC078", "RC078", "RC078", "RC104", "RC104", "TN397", "TN397", "TN397", "TN397", "TN397", "TN412"),
  exp = c("UKG1", "UKG1", "UKH1", "UKH2", "UKG1", "UKH1", "UKH1", "UKH2", "UKH3", "UKH4", "UKH5", "UKG"),
  a_pred = c(300, 400, 200, 350, 200, 600, 125, 200, 300, 500, 125, 125),
  b_pred = c(0.4, 4.1, 0.15, 0.4, 0.6, 0.35, 0.1, 0.1, 0.2, 0.3, 0.6, 0.3)
)



###Pull out dissolved environmental concentrations for each dataset:
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
  separate(KinExp_ID, into=c("cruise", "exp")) %>%
  filter(!exp == "UKG2")






### Create large monte carlo data set by sampling error around each value:

# note: multiplying by sqrt(n-1) (in this case n-1 = 7), 
#  where n = number of points in the calibration curve (in this case 8), converts standard error 
#. to an approximate standard deviation of expected values.

#define number of models
reps <- 1000
set.seed(981)

#generate MC dataset for uptake kinetics experiments
dat.Kin.mc <- dat.MM %>%
  filter(Fragment_mz %in% c(97.0954, 62.8977)) %>%
  group_by(cruise, exp, nM.per.hour.1, se.nM.per.hour.1) %>%
  mutate(mod.val = list(rnorm(reps, mean = nM.per.hour.1,
                              sd = se.nM.per.hour.1 * sqrt(7)))) %>%
  unnest(mod.val) %>%
  mutate(model = paste("model",rep(1:reps),sep="_")) %>%
  ungroup() %>%
  select(cruise, exp, model, treatment_conc, mod.val) %>%
  mutate(mod.val = case_when(mod.val < 0 ~ 0,
                             TRUE ~ mod.val))

#Generate MC dataset for Substrate Concentrations
dat.S.mc <-  exp.enviro.dat %>%
  select(cruise, exp, Mean.Diss.Conc.nM, SD.Diss.Conc.nM) %>%
  group_by(cruise, exp) %>%
  mutate(mod.S.val = list(rnorm(reps, mean = Mean.Diss.Conc.nM,
                                sd = SD.Diss.Conc.nM))) %>%
  unnest(mod.S.val) %>%
  mutate(model = paste("model",rep(1:reps),sep="_")) %>%
  select(cruise, exp, model, mod.S.val) %>%
  mutate(mod.S.val = case_when(mod.S.val < 0 ~ 0,
                               TRUE ~ mod.S.val)) %>%
  filter(!exp == "UKG2")

#Combine MC datasets of kinetics and substrates:
dat.MM.nls <- left_join(dat.Kin.mc, dat.S.mc) %>%
  rename("treatment_nM" = treatment_conc,
         "nM_per_hour" = mod.val,
         "diss_conc_nM" = mod.S.val)  %>%
  unite(col = "exp", cruise:model, remove = FALSE)


#create experiment names lists + prediction matching
mod.names <- dat.MM.nls %>%
  select(model) %>%
  unique()

mod.start.estimates.MM.nls <- model.starting.estiamtes %>%
  cross_join(mod.names) %>%
  select(cruise, exp, model, everything()) %>%
  unite(col = "exp", cruise:model)

dat.exp.names <- mod.start.estimates.MM.nls %>%
  select(exp) %>%
  rename("exp.nls" = exp) %>%
  unique() 





# Run Monte Carlo Model to fit MM equation using NLS ----------------------

#create ouput dataframe
MM.nls.out <- data.frame(
  experiment = as.character(),
  parameter = as.character(),
  Estimate = as.numeric(),
  s_val = as.numeric(),
  sigma = as.numeric(),
  cor = as.numeric(),
  iterations = as.numeric()
)

#fit model (takes a long time)
for (i in dat.exp.names$exp.nls) {
  x <- fit_nls(dat.MM.nls, mod.start.estimates.MM.nls, i)
  MM.nls.out <- rbind(MM.nls.out, x)
}

NLS.output <- MM.nls.out


#Remove outlier models 
NLS.qc <- NLS.output %>%
  separate(experiment, into = c("Cruise", "exp", "text", "model"), remove = FALSE) %>%
  select(experiment, Cruise, exp, model, sigma, cor, iterations) %>%
  unique() %>%
  group_by(Cruise, exp) %>%
  mutate(mean.sigma = mean(sigma),
         sd.sigma = sd(sigma),
         mean.cor = mean(cor),
         sd.cor = sd(cor)) %>%
  mutate(sig.flag.high = sigma > (mean.sigma+(2*sd.sigma)),
         sig.flag.low = sigma < (mean.sigma - (2*sd.sigma)),
         cor.flag.high = cor > (mean.cor + (2*sd.cor)),
         cor.flag.low = cor < (mean.cor - (2*sd.cor)),
         iter.flag = iterations > 50) %>%
  group_by(Cruise, exp, model) %>%
  filter(any(!sig.flag.high == TRUE),
         any(!sig.flag.low == TRUE),
         any(!cor.flag.high == TRUE),
         any(!cor.flag.low == TRUE),
         any(!iter.flag == TRUE))


#Summarize NLS parameters
NLS.sum <- NLS.qc %>%
  left_join(., NLS.output) %>%
  pivot_wider(names_from = parameter, values_from = Estimate) %>%
  separate(experiment, into = c("Cruise", "exp", "text", "model"), remove = FALSE) %>%
  group_by(Cruise, exp) %>%
  reframe(mean_ks = mean(Ks),
          median_ks = median(Ks),
          sd_ks = sd(Ks),
          max_ks = max(Ks),
          min_ks = min(Ks),
          mean_vmax = mean(Vmax),
          sd_vamx = sd(Vmax),
          max_vmax = max(Vmax),
          min_vmax = min(Vmax))





# Calculate Turnover Times and Fluxes using MM Equations ------------------
tt.flux.calcs <- exp.enviro.dat %>%
  rename(Cruise = cruise) %>%
  left_join(., NLS.sum) %>%
  select(Cruise, exp, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, mean_ks, sd_ks, mean_vmax, sd_vamx) %>%
  mutate(v_insitu = (Mean.Diss.Conc.nM*mean_vmax)/(Mean.Diss.Conc.nM+mean_ks),
         v_insitu_sd = v_insitu*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(sd_vamx/mean_vmax)^2+(sd_ks/mean_ks)^2),
         mm_tt = Mean.Diss.Conc.nM/v_insitu,
         mm_tt_sd = mm_tt*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(v_insitu_sd/v_insitu)^2),
         mm_flux_nM_hr = v_insitu,
         mm_flux_sd_nM_hr = v_insitu_sd) %>%
  #  mm_flux = Mean.Diss.Conc.nM/mm_tt,
  #  mm_flux_sd = mm_flux*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(mm_tt_sd/mm_tt)^2)) %>%
  mutate(mm_flux_nM_day = mm_flux_nM_hr*24,
         mm_flux_sd_nM_day = mm_flux_sd_nM_hr*24) %>%
  select(-v_insitu, -v_insitu_sd)



####Organize data and export:
write_csv(tt.flux.calcs, file = "Intermediates/NLS_Kin_Flux_vals.csv")




















































