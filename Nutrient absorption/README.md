# Do fruit secondary metabolites modify nutrient absorption? Insight from fruit bats

# Data information:
The study looks at the impact of different concentration of amides on the nutrient absorption in fruit bats (Carollia perspicilla). The study was conducted by Mariana Gelambi, Susan Whitehead, and Estefania Morales at La Selva Biological Station, Costa Rica during June-July 2021. Bats were fed with a synthetic maintenance diet consisting of mashed ripe banana, soy protein isolate, and supplements. We added piperine (Sigma-Aldrich) to the artificial diet. Four different concentrations of piperine, namely 0.1%, 0.5%, 1.5%, and 2% dry weight were added to the artificial diet, representing the range of natural variation of amides in ripe fruit. We collected fecal samples and quantify carbohydrates and proteins.

# Files: 

## *Code*: nutrient_absorption.R

## *Data*:

(1) raw_chromatogram_data.xlsx: raw data from the peak integration from ChemStation. IT was NOT use in the code.

(2) protein_absorption: information about the experiment, including data, treatment and bat ID.

(3) sugar_absorption.csv: information about the experiment, including data, treatment and bat ID. 

HPLC chromatogram files: 

(4) input_GCalign_sugars.txt: peak area and retention time of each chromatogram. This file was used to align the chromatograms using the package GCalign. The standard format is a tab-delimited text file according to the following layout: (1) The first row contains sample names, the (2) second row column names of the corresponding peak lists. Starting with the third row, peak lists are included for every sample that needs to be incorporated in the data set. Here, a peak list contains data for individual peaks in rows, whereas columns specify variables in the order given in the second row of the text file. Peak lists of individual samples are concatenated horizontally and need to be of the same width (i.e. the same number of columns in consistent order).

(5) output.xlsx and output_trans.csv: Aligned chromatograms and extra info about the samples. The extra info is required to run the models.
