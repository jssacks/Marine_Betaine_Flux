# Marine_Betaine_Flux

The cycling of glycine betaine and homarine in marine microbial communities: quantitative flux measurements and the role of competitive uptake inhibition  


A repository for the data and code for Sacks et al. (in prep). Written in R. 

___

## Dependencies:

In Progress

___ 
# Repo structure:

# Data: 

## Raw_Data:

### Kinetics_Data:
* This folder contains the skyline output .csv files of integrated peaks of targeted liquid chromatography-mass spectrometry data from the uptake kinetics and uptake competition experiments as well as the in situ particulate metabolite measurements of homarine and glycine betaine.

### Distribution_Data:
* This folder contains the skyline output .csv files of integrated peaks of liquid chromatography-mass spectrometry data from the dissolved samples processed using cation-exchange solid phase extraction.


## Meta_Data:

### Data_From_Other_Studies:
* This folder contains uptake rate experimental data from Boysen et al. 2022, 14C incubation primary production estimates for Puget Sound from Newton and VanVoorhis 2002, and a compilation of uptake competition studies from the literature on glycine betaine, DMSP, and choline. 

### Enviro_Meta_Data:
* This folder contains environmental data (particulate carbon, chlorophyll a, temperature, salinity, nutrients, etc.) collected by collaborators or our lab from the research cruises KM1906, TN397, TN412, RC078, and RC104. Most data are also available on Simons CMAP.

### Ingalls_Lab_Data:
* This folder contains data related to sample collection, sample processing, and metabolite quantification. This includes extraction efficiencies from for the cation-exchange solid phase extraction approach, concentrations of standards, and experiment locations. 

___

# Data Processing, Analysis, and Visualization:

## R_Code:

## Folder 1: Meta_Data_Processing

### Script 1.1: 
#### G4_EnviroMeetaData_Interpolation.R
* Sources environmental data (PC, nutrients, chl a, etc.) from TN397 folder in Enviro_Meta_Data, collocates all data based on time, interpolates all data using a Steinman interpolation, and exports a finished, collocated interpolated environmental data csv to the Intermediates folder. 

### Script 1.2:
#### G5_EnviroMetaData_Interpolation.R
* Sources environmental data (PC, nutrients, chl a, etc.) from TN412 folder in Enviro_Meta_Data, collocates all data based on time, interpolates all data using a Steinman interpolation, and exports a finished, collocated interpolated environmental data csv to the Intermediates folder. 

### Script 1.3: 
#### RC078_EnviroMetaData_Curation.R
* Organizes and tidies environmental data from RC078 folder in Enviro_Meta_Data, collocates data based on station and depth, and exports a finished, collocated interpolated environmental data csv to the Intermediates folder.

### Script 1.4:
#### Kinetics_Experiment_MetaData_Curation.R
* Sources data on TN397 (G4), TN412 (G5), KM1906 (G3), RC078 (D1), RC104 (D2), and Newton and VanVoorhis 2002 (PS.PP), matches environmental data against kinetics and uptake experiments and standardizes formatting, and exports a single csv to containing all environmental meta data for each experiment to the Intermediates folder. 

--- 


## Folder 2: Dissolved_Metabolomics_Analysis

### Script 2.1:
#### DissMetab_Functions.R
* Script containing functions used by other scripts to quantify dissolved metabolomics data.  

### Script 2.2: 
#### Raw_Data_Organization.R
* Sources raw data from Raw_Data/Distribution_Data, tidies, organizes, compiles, and exports as a csv to Intermediates.  

### Script 2.3:
#### Dissolved_BMIS.R
* Sources organized raw data from Intermediates folder, runs best-matched internal standard normalization (BMIS) from Boysen et al. (2018), and exports BMISed data to Intermediates folder. 

### Script 2.4:
#### Dissolved_Blk_LOD_Calcs.R
* Sources BMISed data from Intermediates folder, calculates limits of detection (LODs) using laboratory cation-exchange clean blanks (as in Sacks et al. 2022), and exports LODs to Intermediates folder. 

### Script 2.5:
#### Dissolved_QC.R
* Sources organized raw data and LODs from Intermediates folder, runs QC based on changeable thresholds for minimum peak area and LODs, and exports BMISed data to Intermediates folder.

### Script 2.6:
#### Dissolved_RF_RFratios.R
* Sources organized raw data from Intermediates folder and Ingalls lab standards data from Meta_Data folder, calculates response factors (RFs) and response factor ratios (RFratios) as in Boysen et al. (2018), and exports RFs and RFratios as a csv to Intermediates folder. 


### Script 2.7:
#### Dissolved_Quantification.R
* Sources BMISed data, organized raw data, RFs and RFratios, and LODs from Intermediates folder, runs quantification based on RFs and internal standards when applicable, applies an extraction efficiency correction (as in Sacks et al. 2022), and calculates concentrations in samples in nmol/L and nmol C/L and nmol N/L space. 


### Script 2.8:
#### Final_QC_and_Organization.R
* This script sources quantified data and QCed data, removes data not passing QC, and exports data to a csv in Intermediates. 

--- 

## Folder 3: Kinetics_Flux_Analysis.R 

### Script 3.1:
#### Functions.R
* This script contains functions used by other scripts in Kinetics_Flux_Analysis.R

### Script 3.2:
#### Particulate_QE_Quantification.R
* This script sources data of particulate metabolite measurements from Data_Raw/Kinetics_Data, quantifies glycine betaine and homarine using matched internal standards, and exports as a csv to Intermediates. 

### Script 3.3: 
#### Uptake_Rate_Quantification.R
* This script sources data from kinetics experiments from Data_Raw/Kinetics_Data, quantifies the concentration of isotopically labeled glycine betaine and homarine using standard curves, performs blank subtraction on all samples, organizes the experiment incubation time data, and calculates an uptake rate for each sample, and exports this quantified uptake rate data to Intermediates as a csv. 

### Script3.4:
#### UptakeCompetitonAnalysis.R
* This script sources quantified uptake rate data from Intermediates, tests for competitive uptake inhibition of homarine in two uptake competition experiments using Dunett’s test, and exports results to a csv in Intermediates. 

### Script 3.5:
#### NLS_Kinetics_TT_Flux.R
* This script sources quantified uptake rates and dissolved metabolite data from Intermediates, generates Monte Carlo datasets by sampling uncertainty estimates for all values, uses non-linear least-squares to fit Michaelis-Menten (MM) equations, calculates turnover times at in situ concentrations from the MM parameters, calculates mean and SD of MM parameters and turnover times of all non-outlier models, and exports final summarized data to csv in Intermediates. 

### Script 3.6: 
#### WH_Kinetics_TT_Flux.R
* This script sources quantified uptake rates and dissolved metabolite data from Intermediates, generates Monte Carlo datasets by sampling uncertainty estimates for all values, uses Wright-Hobbie linear transformations to determine Michaelis-Menten (MM) parameters, calculates mean and SD of MM parameters and turnover times of all non-outlier models, and exports final summarized data to csv in Intermediates.

### Script 3.7:
#### Kinetics_Flux_Compilation_Context.R
* This script sources kinetics parameters, meta data, particulate and dissolved concentrations from Intermediates, compiles all kinetics and flux results and combines with additional data about those samples, contextualizes flux measurements, calculates Kt/Sn values, and exports all calculated values as csvs to Intermediates. 

### Script 3.8:
#### Statistical_Comparisons.R
* This script sources kinetics parameters, meta data, and Kt/Sn values from Intermediates and runs statistical tests reported in manuscript. 


--- 

Folder 4: Tables

### Script 4.1:
#### Dissolved_Metabolome_Supp_Table.R
* This script organizes dissolved metabolite data and formats them for export. This script produces Supplemental Table 2.

### Script 4.2:
#### Enviro_Metadata_Supplemental_Table.R
* This script organizes environmental meta data and formats them for export. This script produces Supplemental Table 1. 

### Script 4.3:
#### Kin_Flux_Tables.R
* This script organizes final data outputs and formats them for export. This script produces main text Tables 1 and 2 and Supplemental Tables 4-6. 

### Script 4.4:
#### Raw_Data_Tables.R
* This script organizes mass spectrometry peak lists exported from Skyline for dissolved HPLC-MS Q-Exactive orbitrap metabolomics data, particulate HPLC-MS Q-Exactive orbitrap metabolomics data, and single reaction monitoring particulate HPLC-MS TQ-S data and exports to csvs. This script produces Supplemental Tables 9-11. 

___ 

## Folder 5: Figures

## Script 5.1:
### Exp_Location_Maps.R
* This script produces main text Figure 1. 

### Script 5.2:
#### Kinetics_Competition_Experiment_Figure.R
* This script produces main text Figure 2. 

### Script 5.3:
#### Kin_Flux_Conc_Relationships_Figs.R 
* This script produces main text Figure 3. 

### Script 5.4:
#### Blank_Supplemental_Figure_and_Stats.R
* This script produces Supplemental Figure 1 and assess the fit of blank data to a linear model. 

### Script 5.5:
#### NLS_WH_Comparison_Figure.R
* This script produces Supplemental Figure 2. 

### Script 5.6:
#### Particulate_Concentration_Supplemental_Figure.R
* This script produces Supplemental Figure 3. 

### Script 5.7:
#### Environmental_Dissolved_Metabolome_Plot.R 
* This script produces Supplemental Figure 4. 

### Script 5.8:
#### Competition_Network_Lit_Synth_Figure.R 
* This script produces Supplemental Figure 5.

### Script 5.9:
#### All_Kin_Exp_Plot.R 
* This script produces Supplemental Figure 6. 

### Script 5.10:
#### KT_SN_Ratio_PowerLaw.R
* This script produces Supplemental Figure 7. 


--- 

## Output Folders:

## Intermediates:
This folder contains all intermediates files used in workflow. 

## Tables:
This folder contains final tables used in manuscript that are produced by scripts contained in the R_Code/Tables folder. 

## Figures:
This folder contains final figures used in manuscript that are produced by scripts contained in the R_Code/Figures folder. 

___

## Citation

### TBD


![image](https://github.com/user-attachments/assets/2c79f061-37e9-46f6-9ac3-b60c05c8f28a)



![image](https://github.com/user-attachments/assets/67089b18-b436-490b-9969-e7ec920dd006)
