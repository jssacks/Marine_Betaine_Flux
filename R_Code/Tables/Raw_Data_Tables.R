

library(tidyverse)



#define inputs:
tqs.file <- "Intermediates/TQS_dat.csv"
diss.metab.file <- "Intermediates/dissolved_osmo_data_raw.csv"
part.metab.file <- "Intermediates/qe.dat.csv"



###TQS data table organization:
tqs.table <- read_csv(tqs.file) %>%
  select(-date, -type) %>%
  rename("Cruise" = cruise,
         "Experiment" = exp,
         "Treatment" = treatment,
         "Replicate" = rep,
         "Peak_Area" = Area) %>%
  select(SampID, Cruise, Compound, Fragment_mz, Experiment, Treatment, Replicate, Peak_Area)

#export
write_csv(tqs.table, file = "Tables/TQS_RawPeakArea_Table.csv")


###Dissolved raw data file:
diss.dat <- read_csv(diss.metab.file) %>%
  mutate(keep = case_when(str_detect(Rep, "UKH") ~ "Yes",
                          str_detect(Rep, "UKG") ~ "Yes",
                          str_detect(Rep, "UCH") ~ "Yes",
                          str_detect(Rep, "S7_C3_D6") ~ "Yes",
                          str_detect(Rep, "GBT") ~ "Yes",
                          str_detect(Rep, "GBT") ~ "Yes",
                          str_detect(Rep, "Std") ~ "Yes",
                          str_detect(Rep, "Blk") ~ "Yes",
                          str_detect(Rep, "Poo") ~ "Yes",
                          TRUE ~ "No")) %>%
  filter(!keep == "No") %>%
  select(-keep) %>%
  rename("SampID" = Rep)

#export
write_csv(diss.dat, file = "Tables/QE_Diss_RawPeakArea_Table.csv")


###Particulate raw data file 
part.dat <- read_csv(part.metab.file) %>%
  rename("SampID" = Rep)

write_csv(part.dat, file = "Tables/QE_Part_RawPeakArea_Table.csv")







