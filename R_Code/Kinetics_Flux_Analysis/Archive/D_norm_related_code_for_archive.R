



###_____ Section 4: Dnorm Calculations: ___________________________________________________
Kt_S_ratio.dat <- nls.context %>%
  select(Cruise, exp, Compound, mean_ks, sd_ks) %>%
  left_join(., diss.enviro.sum %>% 
              select(Cruise, exp, Compound, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, Rel.Abun, Rel.Abun.SD) %>%
              mutate(Compound = str_replace(Compound, "Glycine betaine", "GBT"))) %>% 
  mutate(kt_s_ratio = mean_ks/Mean.Diss.Conc.nM,
         kt_s_ratio_error = kt_s_ratio*sqrt((SD.Diss.Conc.nM/Mean.Diss.Conc.nM)^2+(sd_ks/mean_ks)^2),
         kt_s_ratio_error_fig_max = kt_s_ratio+kt_s_ratio_error,
         kt_s_ratio_error_fig_min = case_when(kt_s_ratio-kt_s_ratio_error < 0 ~ 0,
                                              kt_s_ratio-kt_s_ratio_error >= 0 ~ kt_s_ratio-kt_s_ratio_error)) %>%
  filter(kt_s_ratio > 1)


Dnorm.dat <- nls.context %>%
  select(Cruise, exp, Compound, mean_ks, sd_ks) %>%
  left_join(., diss.enviro.sum %>% 
              select(Cruise, exp, Compound, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, Rel.Abun, Rel.Abun.SD) %>%
              mutate(Compound = str_replace(Compound, "Glycine betaine", "GBT"))) %>% 
  mutate(D_norm = (mean_ks - Mean.Diss.Conc.nM)/(mean_ks)*100,
         D_norm_2 = mean_ks/Mean.Diss.Conc.nM,
         D_norm_max = (mean_ks+sd_ks - Mean.Diss.Conc.nM-SD.Diss.Conc.nM)/(mean_ks+sd_ks)*100,
         D_norm_min = (mean_ks-sd_ks - Mean.Diss.Conc.nM+SD.Diss.Conc.nM)/(mean_ks-sd_ks)*100) %>%
  filter(D_norm > 0)

####
ggplot(Kt_S_ratio.dat, aes(y = kt_s_ratio, x = Rel.Abun)) +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(x = Rel.Abun, ymin = kt_s_ratio_error_fig_min, ymax = kt_s_ratio_error_fig_max)) +
  geom_errorbarh(aes(y = kt_s_ratio, xmin = Rel.Abun-Rel.Abun.SD, xmax = Rel.Abun+Rel.Abun.SD)) +
  geom_point(aes(color = Compound, size = 2)) +
  scale_y_log10() +
  scale_x_log10()

###
ggplot(Kt_S_ratio.dat, aes(y = kt_s_ratio, x = Rel.Abun)) +
  #  geom_smooth(method = "lm") +
  geom_errorbar(aes(x = Rel.Abun, ymin = kt_s_ratio_error_fig_min, ymax = kt_s_ratio_error_fig_max)) +
  geom_errorbarh(aes(y = kt_s_ratio, xmin = Rel.Abun-Rel.Abun.SD, xmax = Rel.Abun+Rel.Abun.SD)) +
  geom_point(aes(color = Compound, size = 2))# +
#  scale_y_log10() +
# scale_x_log10()


###
ggplot(Kt_S_ratio.dat, aes(y = mean_ks, x = Mean.Diss.Conc.nM)) +
  geom_smooth(method = "lm") +
  # geom_errorbar(aes(x = Rel.Abun, ymin = kt_s_ratio_error_fig_min, ymax = kt_s_ratio_error_fig_max)) +
  # geom_errorbarh(aes(y = kt_s_ratio, xmin = Rel.Abun-Rel.Abun.SD, xmax = Rel.Abun+Rel.Abun.SD)) +
  geom_point(aes(color = Compound, size = 2)) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(.~Compound, scales = "free") +
  geom_abline(intercept = 0, slope = 1)




























