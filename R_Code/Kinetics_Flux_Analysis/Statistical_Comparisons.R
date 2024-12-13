



#load packages:
library(tidyverse)
library(rstatix)


#define inputs:
NLS.file <- "Intermediates/NLS_Kin_Flux_context.csv"
WH.file <- "Intermediates/WH_Kin_Flux_context.csv"
Kt_S_ratio.file <- "Intermediates/Kt_S_ratio_results.csv"
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"

enviro.metab.file <- "Intermediates/Dissolved_Metab_Environment_Summary.csv"
#enviro.metab.file <- "Intermediates/Kinetics_Exp_Enviro_Metabolome_Dat.csv"




###

# Comparisons between GBT and Homarine ------------------------------------


#Perform nonparametric statistics on means of Kt and Vmax comparing GBT and Homarine
flux.kin.dat <- read_csv(NLS.file)

####Perform wilcoxon rank sum test on mean Kt values:
kt.comp.comparison <- flux.kin.dat %>%
  wilcox_test(mean_ks ~ Compound)
kt.comp.comparison

vmax.comp.comparison <- flux.kin.dat %>%
  wilcox_test(mean_vmax ~ Compound)
vmax.comp.comparison

tt.comp.comparison <- flux.kin.dat %>%
  wilcox_test(mm_tt ~ Compound)
tt.comp.comparison

flux.comp.comparison <- flux.kin.dat %>%
  wilcox_test(mm_flux_nM_day ~ Compound)
flux.comp.comparison






# Linear Model Exploring Relationship between Flux and Particulate Concentration --------

#linear model of Flux vs. Part Conc.
log.flux.pconc <- flux.kin.dat %>%
  select(Cruise, exp, Compound, mm_flux_nM_day, mean.part.conc.nM) %>%
  mutate(log10_flux = log10(mm_flux_nM_day),
         log10_pconc = log10(mean.part.conc.nM))


flux.pconc.model <- lm(log10_flux~log10_pconc, data = log.flux.pconc)
summary(flux.pconc.model)


##Just homarine:
hom.flux.pconc.dat <- log.flux.pconc %>%
  filter(Compound == "Homarine")

hom.flux.pconc.model <- lm(log10_flux~log10_pconc, data = hom.flux.pconc.dat)
summary(hom.flux.pconc.model)

##Just GBT:
gbt.flux.pconc.dat <- log.flux.pconc %>%
  filter(Compound == "GBT")

gbt.flux.pconc.model <- lm(log10_flux~log10_pconc, data = gbt.flux.pconc.dat)
summary(gbt.flux.pconc.model)







# Regression of Kt vs. dissolved concentration  ---------------------------

###############
log.kt.dconc.dat <- flux.kin.dat %>%
  mutate(log10_ks = log10(mean_ks),
         log10_dconc = log10(Mean.Diss.Conc.nM))


###Overall model 
kt.dconc.lm <- lm(log10_dconc ~ log10_ks, data = log.kt.dconc.dat)
summary(kt.dconc.lm)


###Just homarine
hom.kt.dconc.dat <- log.kt.dconc.dat %>%
  filter(Compound == "Homarine")

hom.kt.dconc.lm <- lm(log10_dconc ~ log10_ks, data = hom.kt.dconc.dat)
summary(hom.kt.dconc.lm)


###Just GBT
gbt.kt.dconc.dat <- log.kt.dconc.dat %>%
  filter(Compound == "GBT")

gbt.kt.dconc.lm <- lm(log10_dconc ~ log10_ks, data = gbt.kt.dconc.dat)
summary(gbt.kt.dconc.lm)







# Kt_Sn_Ratio Comparison with % of Dissolved Pool -------------------------

#organize and compile dissolved data and kt_s_ratio data
kt_s_ratio.dat <- read_csv(Kt_S_ratio.file) %>%
  select(Cruise, exp, Compound, mean_ks, sd_ks, kt_s_ratio, kt_s_ratio_error)

dissolved.metab.dat <- read_csv(enviro.metab.file) %>%
  select(Cruise, exp, Compound, Mean.Diss.Conc.nM, SD.Diss.Conc.nM) %>%
  filter(!Compound == "Trimethylamine N-oxide") %>%
  mutate(Compound = str_replace(Compound, "Glycine betaine", "GBT")) %>% 
  group_by(Cruise, exp) %>%
  mutate(Tot.Sulf.Bet.Conc.nM = sum(Mean.Diss.Conc.nM),
         Tot.Sulf.Bet.Conc.nM.SD = sum(SD.Diss.Conc.nM),
         Perc.tot.pool = Mean.Diss.Conc.nM/Tot.Sulf.Bet.Conc.nM*100,
         Perc.tot.pool.SD =  Perc.tot.pool*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(Tot.Sulf.Bet.Conc.nM.SD/Tot.Sulf.Bet.Conc.nM)^2))

kt_s_ratio.test.dat <- left_join(kt_s_ratio.dat, dissolved.metab.dat) %>%
  mutate(log10.ktsratio = log10(kt_s_ratio),
         log10.perctotpool = log10(Perc.tot.pool))


##run linear model on log10 transformed data:
###############

###Overall model 
ktsratio.lm <- lm(log10.ktsratio ~ log10.perctotpool, data = kt_s_ratio.test.dat)
summary(ktsratio.lm)


###Just homarine
hom.kt_s_ratio.test.dat <- kt_s_ratio.test.dat %>%
  filter(Compound == "Homarine")

hom.ktsratio.lm <- lm(log10.ktsratio ~ log10.perctotpool, data = hom.kt_s_ratio.test.dat)
summary(hom.ktsratio.lm)


###Just GBT
gbt.kt_s_ratio.test.dat <- kt_s_ratio.test.dat %>%
  filter(Compound == "GBT")

gbt.ktsratio.lm <- lm(log10.ktsratio ~ log10.perctotpool, data = gbt.kt_s_ratio.test.dat)
summary(gbt.ktsratio.lm)


#Comparison of total pool with Kt
ggplot(kt_s_ratio.test.dat, aes(x = Tot.Sulf.Bet.Conc.nM, y = mean_ks)) +
  geom_smooth(method = "lm", alpha = 0.3) +
  geom_point(size = 2, aes(color = Compound)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_y_log10() +
  scale_x_log10() +
  theme_test() +
  xlab("Total Betaine + Sulfonium Pool (nM)") +
  ylab("Kt (nM)")


###Overall model 
totpool.kt.lm <- lm(log10(mean_ks) ~ log10(Tot.Sulf.Bet.Conc.nM), data = kt_s_ratio.test.dat)
summary(totpool.kt.lm)








# Comparison of NLS and WH values -----------------------------------------

#load in data and rename variables based on approach used:

#NLS
nls.dat <- read_csv(NLS.file) %>%
  mutate(nls.ktsn = mean_ks+Mean.Diss.Conc.nM,
         nls.ktsn.sd = sd_ks+SD.Diss.Conc.nM,
         nls.vmax = mean_vmax,
         nls.vmax.sd = sd_vamx,
         nls.tt = mm_tt,
         nls.tt.sd = mm_tt_sd) %>%
    select(Cruise, exp, Compound, nls.ktsn, nls.ktsn.sd, nls.vmax, nls.vmax.sd, nls.tt, nls.tt.sd)
  
#Wright-Hobbie
wh.dat <- read_csv(WH.file) %>%
  mutate(wh.ktsn = wh_mean_kt_sn,
         wh.ktsn.sd = wh_sd_kt_sn,
         wh.vmax = wh_mean_vmax,
         wh.vmax.sd = wh_sd_vmax,
         wh.tt = wh_mean_tt,
         wh.tt.sd = wh_sd_tt) %>%
  select(Cruise, exp, Compound, wh.ktsn, wh.ktsn.sd, wh.vmax, wh.vmax.sd, wh.tt, wh.tt.sd)


##combine datasets: 
nls.wh.dat <- left_join(nls.dat, wh.dat)



#Kt+Sn__________

#visualize: 
ggplot(nls.wh.dat, aes(x = wh.ktsn, y = nls.ktsn)) + 
  geom_smooth(method = "lm", alpha = 0.3) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

#linear model:
ktsn.lm <- lm(nls.ktsn~wh.ktsn, data = nls.wh.dat)
summary(ktsn.lm)



### Vmax_________

#visualize:
ggplot(nls.wh.dat, aes(x = wh.vmax, y = nls.vmax)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) 

#linear model:
vmax.lm <- lm(nls.vmax~wh.vmax, data = nls.wh.dat)
summary(vmax.lm)



#### TT_________

#visualize
ggplot(nls.wh.dat, aes(x = wh.tt, y = nls.tt)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

#linear model:
tt.lm <- lm(nls.tt~wh.tt, data = nls.wh.dat)
summary(tt.lm)

































