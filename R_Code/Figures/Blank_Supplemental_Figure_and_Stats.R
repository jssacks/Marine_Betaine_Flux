


#packages:
library(tidyverse)
library(ggsci)


#inputs
blk.file <- "Intermediates/Kinetics_Blanks_Final.csv"



##load in data and tidy dataframe, 
dat.blk <- read_csv(blk.file) %>%
  mutate(Compound = case_when(Fragment_mz == 45.1882 ~ "GBT",
                              Fragment_mz == 62.8977 ~ "GBT",
                              TRUE ~ "Homarine")) %>%
  filter(Fragment_mz %in% c(62.8977, 78.0399)) 


##Run linear models to estimate significance + R2 values
dat.hom <- dat.blk %>%
  filter(Compound == "Homarine")

dat.gbt <- dat.blk %>%
  filter(Compound == "GBT")


##homarine
hom.lm <- lm(nM_in_vial~treatment, data = dat.hom)
summary(hom.lm)

##GBT
gbt.lm <- lm(nM_in_vial~treatment, data = dat.gbt)
summary(gbt.lm)



#make linearly scaled plot
plot.lin <- ggplot(dat.blk, aes(x = treatment, y = nM_in_vial)) +
  geom_smooth(method = "lm", color = "black", alpha = 0.3) +
  geom_point(shape = 21, size = 3, stroke = 0.2, aes(fill = cruise)) +
  facet_wrap(.~Compound, scales = "free") +
  xlab(expression("Treatment Concentration"~(nmol~L^-1))) +
  ylab(expression("Concentration in Vial"~(nmol~L^-1))) +
  theme_bw()
plot.lin

#make log scaled plot
plot.log <- ggplot(dat.blk, aes(x = treatment, y = nM_in_vial)) +
  geom_smooth(method = "lm", color = "black", alpha = 0.3) +
  geom_point(shape = 21, size = 3, stroke = 0.2, aes(fill = cruise)) +
  facet_wrap(.~Compound, scales = "free") +
  scale_y_log10() +
  scale_x_log10() +
  xlab(expression("Treatment Concentration"~(nmol~L^-1))) +
  ylab(expression("Concentration in Vial"~(nmol~L^-1))) +
  theme_bw()
plot.log

blk.plot <- ggarrange(plot.lin,
                      NA,
                      plot.log,
                      nrow = 3, ncol = 1,
                      heights = c(0.48, 0.04, 0.48),
                      labels = c("A", NA, "B"),
                      common.legend = TRUE,
                      legend = "bottom")
blk.plot

ggsave(blk.plot,  filename = "Figures/SupFig1_Blanks.png",
       height = 5, width = 5, dpi = 300, units = "in", scale = 1.3, bg = "white")




















































