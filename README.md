# Interactions between nutrients and fruit secondary metabolites shape bat foraging behavior and nutrient absorption

This repository contains scripts, data files, and outputs related to the research paper titled "Interactions between Nutrients and Fruit Secondary Metabolites: Bat Foraging Behavior and Nutrient Absorption." The study looks at the impact of different concentration of piperine on the foraging behavior and nutrient absorption in fruit bats (*Carollia perspicilla*). The study was conducted by Mariana Gelambi, Susan Whitehead, and Estefania Morales at La Selva Biological Station, Costa Rica during June-July 2021. Bats were fed with a synthetic maintenance diet consisting of mashed ripe banana, soy protein isolate, and supplements. We added piperine (Sigma-Aldrich) to the artificial diet. Four different concentrations of piperine, namely 0.1%, 0.5%, 1.5%, and 2% dry weight were added to the artificial diet, representing the range of natural variation of amides in ripe fruit. We collected fecal samples and quantify carbohydrates and proteins using HPLC and spectrophotometric assays, respectevly. 

## Objective 1. The relative role of nutrients and defensive metabolites in bat preference

### Scripts

#### 1. `bat_preference.R`

This script analyzes preference experiments using Generalized Linear Mixed Models (GLMMs).

### Data Files

#### 1. `preference.csv`

- `Date`: Date in mm/dd/yyyy format.
- `Bat`: Bat identification.
- `BatID`: A unique number for each individual bat.
- `Treatment`: Treatment categories:
  - `ha`: High amides
- `hn`: High nutrients
- `hn+hamides`: High amides + high nutrients
- `la`: Low amides
- `ln`: Low nutrients
- `ln+lamides`: Low amides + low nutrients

## Objective 2. The effect of piperine on sugar and protein absorption

### Scripts

#### 1. `nutrient_absorption.R`

This script processes raw chromatograms, analyzes protein and sugar absorption using GLMMs and NMDS.

### Data Files

#### 1. `raw_chromatogram_data.xlsx`

Raw data from peak integration from ChemStation.

#### 2. `protein_absorption.csv`

Experiment information, including data, treatment, and bat ID for protein absorption analysis.

#### 3. `sugar_absorption.csv`

Experiment information, including data, treatment, and bat ID for sugar absorption analysis.

#### 4. HPLC Chromatogram Files:

- `input_GCalign_sugars.txt`: Peak area and retention time of each chromatogram for chromatogram alignment using the GCalign package (tab-delimited text file).
- `output.xlsx` and `output_trans.csv`: Aligned chromatograms and metadata.

## Plots and Tables

The 'Plots and Tables' folder contains various output files generated from the analyses conducted in the scripts. These files visually represent the results and insights obtained from the data exploration and statistical modeling.
