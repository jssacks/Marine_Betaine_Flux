




##load packages
library(tidyverse)

##Functions:
source("R_Code/Dissolved_Metabolomics_Analysis/DissMetab_Functions.R")


#
#_________________Dissolved_____________
#G3 Dissolved
g3d.file <- "Raw_Data/Distribution_Data/Dissolved/G3_dissolved_HILICPos_Oct24.csv"

#G4 Dissolved
g4d.file.0 <-  "Raw_Data/Distribution_Data/Dissolved/G4_dissolved_HILICPos_Oct24_pooled.csv"
g4d.file.1 <-  "Raw_Data/Distribution_Data/Dissolved/G4_dissolved_HILICPos_Oct24_file1.csv"
g4d.file.2 <-  "Raw_Data/Distribution_Data/Dissolved/G4_dissolved_HILICPos_Oct24_file2.csv"

#D1 Dissolved
d1d.file <- "Raw_Data/Distribution_Data/Dissolved/D1_dissolved_HILICPos_Oct24.csv"

#Kinetics Experiments Dissolved
ked.file <-  "Raw_Data/Distribution_Data/Dissolved/KinExp_dissolved_HILICPos_Oct24.csv"



# Organize Dissolved Data ---------------------------------------------

#Organize G3-dissolved
g3d <- sky_read(g3d.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "KM1906")

#Organize G4-dissolved
g4d.0 <- sky_read(g4d.file.0) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "TN397")

g4d.1 <- sky_read(g4d.file.1) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "TN397")  %>%
  filter(!str_detect(Rep, "_Std_"))

g4d.2 <- sky_read(g4d.file.2) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "TN397")%>%
  filter(!str_detect(Rep, "_Std_"))

g4d <- rbind(g4d.0, g4d.1, g4d.2)


###Organize D1 Dissolved
d1d <- sky_read(d1d.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "RC078")


#Organize Kinetics Experiment Dissolved
ked <- sky_read(ked.file) %>%
  select(Rep, Compound, Area) %>%
  mutate(Cruise = "KinExp") 

####Pull together all dissolved data:
all.d.dat <- rbind(g3d, g4d, ked, d1d)

###Grab Just IS:
is.d.dat <- all.d.dat %>%
  filter(str_detect(Compound, ", ")) %>%
  filter(!Compound %in% c("Cys-Gly, oxidized"))


###Grab Just osmolytes: 
osmo.d.dat <- all.d.dat %>%
  filter(Compound %in% 
           c("beta-Alaninebetaine", "Glycine betaine", "Proline betaine", "Homarine",
             "Trigonelline", "Betonicine", "Dimethylsulfonioacetate", "Dimethylsulfoniopropionate",
             "Gonyol", "Trimethylamine N-oxide"))

#make sample list:
d.smp.list <- osmo.d.dat %>%
  select(Rep, Cruise) %>%
  unique() %>%
  mutate(Injec_vol = case_when(str_detect(Rep, "Half") ~ 1,
                             TRUE ~ 2))

####Write cleaned up IS and betaine data to .csvs 
write_csv(is.d.dat, file = "Intermediates/dissolved_IS_data_raw.csv")

write_csv(osmo.d.dat, file = "Intermediates/dissolved_osmo_data_raw.csv")

write_csv(d.smp.list, file = "Intermediates/dissolved_smp_list.csv")














