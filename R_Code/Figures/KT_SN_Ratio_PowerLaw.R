



#packages:
library(tidyverse)
library(ggpubr)
library(cowplot)



#define inputs:
flux.kin.file <- "Intermediates/NLS_Kin_Flux_vals.csv"
context.file <- "Intermediates/NLS_Kin_Flux_context.csv"
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"
enviro.metab.file <- "Intermediates/Dissolved_Metab_Environment_Summary.csv"
Kt_S_ratio.file <- "Intermediates/Kt_S_ratio_results.csv"

#region.pal <- c("#00887d", "#014d64", "#01a2d9", "#6794a7")



#load in and combine datasets:

#load in data:

#region classifications
region.dat <- read_csv(meta.data.file) %>%
  select(KinExp_ID, Region) %>%
  separate(KinExp_ID, into = c("Cruise", "exp"), remove = FALSE) %>%
  mutate(Region = case_when(Region == "Gyre" ~ "NPSG",
                            Region == "Coastal" ~ "PS",
                            Region == "Equatorial" ~ "Eq",
                            TRUE ~ Region)) %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("NPSG", "Eq", "PS", "NPTZ"))) 

##Remove TMAO and calculate total betaine + sulfonium concentration:
summed.comp.dat <- read_csv(enviro.metab.file) %>%
  filter(!Compound == "TMAO") %>%
  group_by(Cruise, exp) %>%
  reframe(metab.tot.nM = sum(Mean.Diss.Conc.nM),
          metab.tot.nM.uncert = sqrt(sum(SD.Diss.Conc.nM)^2))

#flux and kinetics data
flux.kin.dat <- read_csv(flux.kin.file) %>%
  select(Cruise, exp, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, mean_ks, sd_ks, mm_flux_nM_day, mm_flux_sd_nM_day)

#particulate concentration data
context.dat <- read_csv(context.file) %>%
  select(Cruise, exp, Compound, mean.part.conc.nM, sd.part.conc.nM)

#Kt_S_ratio data
kt_s_ratio.dat <- read_csv(Kt_S_ratio.file)


##Combine all datasets
all.dat <- flux.kin.dat %>%
  left_join(., context.dat) %>%
  left_join(., summed.comp.dat) %>%
  left_join(., region.dat) %>%
  left_join(., kt_s_ratio.dat) %>%
  mutate(perc_diss_pool = (Mean.Diss.Conc.nM/metab.tot.nM)*100,
         error_perc_diss_pool = (Mean.Diss.Conc.nM*sqrt(((metab.tot.nM.uncert-SD.Diss.Conc.nM)/metab.tot.nM)^2 + (SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2)/metab.tot.nM*100),
         error_perc_diss_pool_fig_max = perc_diss_pool+error_perc_diss_pool,
         error_perc_diss_pool_fig_min = case_when(perc_diss_pool-error_perc_diss_pool < 0.5 ~ 0.5,
                                                  perc_diss_pool-error_perc_diss_pool >= 0.5 ~ perc_diss_pool-error_perc_diss_pool))%>%
  mutate(Compound = as.factor(Compound)) %>%
  filter(!is.na(kt_s_ratio))



##Plot Kt/Sn vs. % of Dissolved Pool:

######Relative Kt difference plotted against % of total betaine + sulfonium pool


plot.1 <- ggplot(all.dat, aes(x = perc_diss_pool, y = kt_s_ratio)) + 
  geom_point(aes(shape = Compound, fill = Region), size = 3.5) +
  scale_shape_manual(values = c(21, 24), drop = FALSE) +
  scale_fill_manual(values = c("NPTZ" =  "#3C5488FF",
                               "NPSG" =  "#E64B35FF",
                               "Eq" = "#4DBBD5FF",
                               "PS" = "#00A087FF"), drop = FALSE) +
  theme_test()  +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  xlab("Percentage of Dissolved Pool (%)") +
  ylab(expression(K[t]/S[n])) 
plot.1

ggsave(plot.1, filename = "Figures/SupFig7_KT_SN_Ratio_PowerLaw.png",
       scale = 1.1, units = "in", dpi = 600,
       width = 6, height = 4.5, bg = "white")











































