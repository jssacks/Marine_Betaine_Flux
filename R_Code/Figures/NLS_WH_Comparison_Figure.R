



#load packages:
library(tidyverse)
library(ggpubr)


#define inputs:
NLS.file <- "Intermediates/NLS_Kin_Flux_context.csv"
WH.file <- "Intermediates/WH_Kin_Flux_context.csv"
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"



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

##load in metadata:
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






##combine datasets: 
nls.wh.dat <- left_join(nls.dat, wh.dat) %>%
  left_join(., region.dat)




###Visualize Kt+Sn data
ktsn.plot <- ggplot(nls.wh.dat, aes(x = wh.ktsn, y = nls.ktsn)) +
  geom_smooth(alpha = 0.2, method = "lm", color = "gray20") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = nls.ktsn-nls.ktsn.sd, ymax = nls.ktsn+nls.ktsn.sd), width = 0.03) +
  geom_errorbarh(aes(xmin = wh.ktsn-wh.ktsn.sd, xmax = wh.ktsn+wh.ktsn.sd), height = 0.03) +
  geom_point(aes(shape = Compound, fill = Region), size = 3.5) +
  scale_shape_manual(values = c(21, 24), drop = FALSE) +
  scale_fill_manual(values = c("NPTZ" =  "#3C5488FF",
                               "NPSG" =  "#E64B35FF",
                               "Eq" = "#4DBBD5FF",
                               "PS" = "#00A087FF"), drop = FALSE) +
  theme_test()  +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  ylab(expression(NLS~K[t]+S[n]~(nmol~L^-1))) +
  xlab(expression(WH~K[t]+S[n]~(nmol~L^-1))) +
  annotate("text", x = 80, y = 500, 
           label = "p < 0.0001",
           size = 4) +
  annotate("text", x = 80, y = 430, 
           label = expression(R^2~"="~0.91),
           size = 4) 
ktsn.plot



###Visualize Vmax data
vmax.plot <- ggplot(nls.wh.dat, aes(x = wh.vmax, y = nls.vmax)) +
  geom_smooth(alpha = 0.2, method = "lm", color = "gray20") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = nls.vmax-nls.vmax.sd, ymax = nls.vmax+nls.vmax.sd), width = 0.03) +
  geom_errorbarh(aes(xmin = wh.vmax-wh.vmax.sd, xmax = wh.vmax+wh.vmax.sd), height = 0.03) +
  geom_point(aes(shape = Compound, fill = Region), size = 3.5) +
  scale_shape_manual(values = c(21, 24), drop = FALSE) +
  scale_fill_manual(values = c("NPTZ" =  "#3C5488FF",
                               "NPSG" =  "#E64B35FF",
                               "Eq" = "#4DBBD5FF",
                               "PS" = "#00A087FF"), drop = FALSE) +
  theme_test()  +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  ylab(expression(NLS~V[max]~(nmol~L^-1~h^-1))) +
  xlab(expression(WH~V[max]~(nmol~L^-1~h^-1))) +
  annotate("text", x = 1, y = 3, 
           label = "p < 0.0001",
           size = 4) +
  annotate("text", x = 1, y = 2.65, 
           label = expression(R^2~"="~0.99),
           size = 4) 
vmax.plot



###Visualize Turnover time data
tt.plot <- ggplot(nls.wh.dat, aes(x = wh.tt, y = nls.tt)) +
  geom_smooth(alpha = 0.2, method = "lm", color = "gray20") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = nls.tt-nls.tt.sd, ymax = nls.tt+nls.tt.sd), width = 0.03) +
  geom_errorbarh(aes(xmin = wh.tt-wh.tt.sd, xmax = wh.tt+wh.tt.sd), height = 0.03) +
  geom_point(aes(shape = Compound, fill = Region), size = 3.5) +
  scale_shape_manual(values = c(21, 24), drop = FALSE) +
  scale_fill_manual(values = c("NPTZ" =  "#3C5488FF",
                               "NPSG" =  "#E64B35FF",
                               "Eq" = "#4DBBD5FF",
                               "PS" = "#00A087FF"), drop = FALSE) +
  theme_test()  +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  ylab("NLS TT (h)") +
  xlab("WH TT (h)") +
  annotate("text", x = 300, y = 1750, 
           label = "p < 0.0001",
           size = 4) +
  annotate("text", x = 300, y = 1500, 
           label = expression(R^2~"="~0.92),
           size = 4) 
tt.plot



###combine all three plots:

#get legend
legend <- get_legend(ktsn.plot)

#remove legend from all plots
ktsn.plot.2 <- ktsn.plot +
  theme(legend.position = "none")

vmax.plot.2 <- vmax.plot +
  theme(legend.position = "none")

tt.plot.2 <- tt.plot +
  theme(legend.position = "none")

plot.all <- ggarrange(NA, NA, NA, NA, NA,
                      NA, ktsn.plot.2, NA, vmax.plot.2, NA,
                      NA, NA, NA, NA, NA,
                      NA, tt.plot.2, NA, legend, NA,
                      NA, NA, NA, NA, NA,
                      nrow = 5, ncol = 5,
                      widths = c(0.02, 0.43, 0.02, 0.43, 0.02),
                      heights = c(0.02, 0.43, 0.02, 0.43, 0.02),
                      labels = c(NA, NA, NA, NA, NA,
                                 NA, "A", NA, "B", NA,
                                 NA, NA, NA, NA, NA,
                                 NA, "C", NA, NA, NA,
                                 NA, NA, NA, NA, NA))

plot.all
ggsave(plot.all, filename = "Figures/SupFig6_WH.NLS.Comparison.png",
       scale = 1.2, units = "in", dpi = 600,
       width = 7, height = 6, bg = "white")






























