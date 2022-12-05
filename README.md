# Interactions between nutrients and fruit secondary metabolites shape bat foraging behavior and nutrient absorption

# Data information:
The study looks at the impact of different concentration of piperine on the foraging behavior and nutrient absorption in fruit bats (*Carollia perspicilla*). The study was conducted by Mariana Gelambi, Susan Whitehead, and Estefania Morales at La Selva Biological Station, Costa Rica during June-July 2021. Bats were fed with a synthetic maintenance diet consisting of mashed ripe banana, soy protein isolate, and supplements. We added piperine (Sigma-Aldrich) to the artificial diet. Four different concentrations of piperine, namely 0.1%, 0.5%, 1.5%, and 2% dry weight were added to the artificial diet, representing the range of natural variation of amides in ripe fruit. We collected fecal samples and quantify carbohydrates and proteins using HPLC and spectrophotometric assays, respectevly. 

# Preference experiments

## Scripts

### 1. bat_preference.R
This script was used to analyze preference experiments using GLMMs. 

## Data files

### 1. preference.csv 
Date in mm/dd/yyyy format.
Bat, batID, one unique number for each individual
Treatment,
ha = high amides
hn = high nutrients
hn+hamides = high amides + high nutrients
la = low amides
ln = low nutrients
ln+lamides = low amides + low nutrients

# Nutrient absorption 

## Scripts  

### 1. nutrient_absorption.R
This script was used to process the raw chromatograms and to analyze the absorption of proteins and sugars using GLMMs and NMDS. 

## Data files

### 1. raw_chromatogram_data.xlsx
Raw data from the peak integration from ChemStation.

### 2. protein_absorption.csv
Information about the experiment, including data, treatment and bat ID.

### 3. sugar_absorption.csv
Information about the experiment, including data, treatment and bat ID. 

### 4. HPLC chromatogram files: 
input_GCalign_sugars.txt: peak area and retention time of each chromatogram. This file was used to align the chromatograms using the package GCalign. The standard format is a tab-delimited text file.
output.xlsx and output_trans.csv: Aligned chromatograms and extra info about the samples. The extra info is required to run the models.

# Plots and Table folders
Contain different outputs from the analyses conducted in the scripts. 

