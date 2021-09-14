# Malpighiaceae_supplementary

This repository contains supplementary files related to the manuscript entitled:
"Untargeted metabolomics sheds light on the characterization and diversity of major classes of secondary metabolites in Malpighiaceae".

# Citation
Helena Mannochio-Russo, Rafael Felipe de Almeida, Wilhan Donizete Gonçalves Nunes, Paula Carolina Pires Bueno, Andres Mauricio Caraballo-Rodriguez, Anelize Bauermeister, Pieter C. Dorrestein, Vanderlan S. Bolzani.
Manuscript in preparation.

# Folders

## Qemistree
This folder contains .qzv files for tree visualization, and the Jupyter Notebook for running Qemistree in Qiime2 plugin (more information in https://github.com/biocore/q2-qemistree).

##### Empress visualization
The .zip contains two .qzv files. These .qzv files are the results obtained from the chemical hierarchy analysis (Qemistree), which can be interactively explored in EMPRESS.
To explore the tree results, visit https://view.qiime2.org/ and drag and drop a .qzv file in the marked box. The tree layout, colors, and metadata information
can all be interactively visualized with the menu on the right. For more information about EMPRESS visualization, please visit https://github.com/biocore/empress.

## Feature classification counts
Contains the CANOPUS classifications for each feature in the positive and negative ionization mode and a Jupyter Notebook used to build a feature table with classification counts (cutoff 1000).

## Statistical analysis
Contains an R script used to determine the statistically significant classes in POS and NEG, as well as two .rds files (POS and NEG) containing the feature table, the metadata information, and the CANOPUS classifications obtained.

## Heatmaps
Contains a CANOPUS classlist containing the statistically relevant classes, a Malpighiaceae phylogeny file, a feature table (cutoff1000) with the statistically significant classes, and a jupyter notebook for the heatmaps construction.

## Formating tables for MLE - presence/absence

