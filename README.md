# Genome-Analysis

Welcome to the GENOME-ANALYSES repository!

This repository contains the code for MDPI article titled ["Compression-complexity measures for analysis and classification of coronaviruses"](https://doi.org/10.3390/e25010081) published in Entropy 2022 by Naga Venkata Trinath Sai Munagala 1†, Prem Kumar Amanchi 1†, Karthi Balasubramanian 1,* Athira Panicker and Nithin Nagaraj. † These authors contributed equally to this work. * Author to whom correspondence should be addressed.

The study focuses on utilizing compression-complexity-based distance measures to analyze genomic sequences, generate phylogenetic trees, and train machine learning models using complexity values as features.

## Setup
These instructions will guide you in setting up the GENOME-ANALYSES project on your local machine.

Place the cloned repo GENOME-ANALYSES in your MATLAB working directory and add it to the path (right-click on the GENOME-ANALYSES directory and select "Add to Path > Folders and Subfolders").
Change the MATLAB current working directory to ~/MATLAB/GENOME-ANALYSES/ (simply open this directory in MATLAB).
Running the Software
## Classification
To run the classification program:

Open the workspace directory GENOME-ANALYSES.
Click the "Run" button in the MATLAB editor.
Alternatively, you can run it from the command window by typing Classification_main (without quotes) and pressing ENTER.

## Phylogeny Tree Generation
To generate the phylogeny tree:

Open the workspace directory GENOME-ANALYSES.
Click the "Run" button in the MATLAB editor.
Alternatively, you can run it from the command window by typing phylogeny_main (without quotes) and pressing ENTER.

## Common Instructions for Classification
Classification_main.m is the main program that you need to run for obtaining classification results.
The complete list of provided datasets used in the paper can be found in the "DataBase" section of the Readme file.
The default dataset is 'covid'. You can change the dataset by modifying line number 8 in Classification_main.m.
If you want to use your own datasets, make sure they consist of at least 2 clusters or subfolders in the DataBase folder. For .fasta files, create subfolders representing different clusters in the directory GENOME-ANALYSES/DataBase/TestDataSet and place the respective sequences (.fasta files) in the subfolders. For NCBI accession numbers, create a .txt file with comma-separated NCBI accession numbers for each cluster and place all .txt files (one per cluster) in the directory GENOME-ANALYSES/DataBase/TestDataSet.
You can choose the feature to be used for classification from the following options:
LZ (Lempel-Ziv complexity)
ETC (Effort-to-compress complexity)
ETC+LZ (Effort-to-compress complexity and Lempel-Ziv complexity)
## Common Instructions for Phylogenetic Tree
phylogeny_main.m is the main program that you need to run for generating the phylogenetic tree.
The complete list of provided datasets used in the paper can be found in the "DataBase" section of the Readme file.
The default dataset is 'Mammals'. You can change the dataset by modifying line number 9 in phylogeny_main.m.
