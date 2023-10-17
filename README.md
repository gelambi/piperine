# Interactions between nutrients and fruit secondary metabolites shape bat foraging behavior and nutrient absorption

This repository contains scripts, data files, and outputs related to the research paper titled "Interactions between Nutrients and Fruit Secondary Metabolites: Bat Foraging Behavior and Nutrient Absorption." The study looks at the impact of different concentration of piperine on the foraging behavior and nutrient absorption in fruit bats (*Carollia perspicilla*). The study was conducted by Mariana Gelambi, Susan Whitehead, and Estefania Morales at La Selva Biological Station, Costa Rica during June-July 2021. Bats were fed with a synthetic maintenance diet consisting of mashed ripe banana, soy protein isolate, and supplements. We added piperine (Sigma-Aldrich) to the artificial diet. Four different concentrations of piperine, namely 0.1%, 0.5%, 1.5%, and 2% dry weight were added to the artificial diet, representing the range of natural variation of amides in ripe fruit. We collected fecal samples and quantify carbohydrates and proteins using HPLC and spectrophotometric assays, respectively. 

## Objective 1. The relative role of nutrients and defensive metabolites in bat preference

### Scripts

#### 1. `script1_objective1.R`

This script addresses Objective 1.1, 1.2, 1.3, and 1.4. For the first three objectives, we conducted a comparison of the amount of food consumed on each Petri dish per bat using paired t-tests, one test per trial (1.1, 1.2, and 1.3). We used the t.test() function, where the null hypothesis was set as the mean difference being equal to 0, indicating no preference between the two options being compared. Then, for Objective 1.4, our focus shifted towards estimating the amount of nutrients (mg) and piperine consumed per bat in each trial. To achieve this, we fitted two separate linear models for nutrients and piperine consumption across the three trials. We used the lm() function to build these models, allowing us to estimate the differences in nutrient and piperine intake between the trials. 

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

This script addresses Objective 2, which focuses on assessing the effects of piperine on the excretion of (2.1) individual and total sugars and (2.2) total proteins.

In the case of sugar excretion, we aligned the peaks based on the retention time using the R package GCalignR (Ottensmann et al. 2018). We calculated the total amount of soluble sugars by adding the areas of all the 12 individual sugar peaks found in the sample and expressed the concentration as glucose equivalents based on a standard curve for glucose. Then, we fitted generalized linear mixed models (GLMMs) to estimate the differences in sugar excretion (individual and total sugars) between the control group and the four concentrations of piperine tested (0.1%, 0.5%, 1%, and 2%). The response variable in these models was the proportion of total sugars per unit fecal mass excreted by bats, while the predictor variable was the concentration of piperine. To account for repeated measurements on the same bats and the use of a fresh artificial diet each evening, we included random effects for bat identity and trial date. To contrast each piperine concentration with the control, we utilized the emmeans() function from the emmeans package.

Similar to Objective 2.1, for protein excretion, we fitted separate GLMMs to estimate the differences in protein excretion between the control group and the four concentrations of piperine. The response variable was the proportion of total proteins per unit fecal mass excreted, and the predictor variable was the piperine concentration. Random effects for bat identity and trial date were also included to account for repeated measurements and dietary variations. The emmeans() function was used to obtain estimated marginal means for each treatment level, enabling a comparison between each piperine concentration and the control.

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

## Figures

The 'Figure' folder contains four figures generated from the analyses conducted in the scripts.