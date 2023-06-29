# Genome-Analysis

Welcome to the GENOME-ANALYSES repository! This project is based on the research paper titled "Compression-complexity measures for analysis and classification of coronaviruses." (DOI: 10.3390/e25010081) The study focuses on utilizing compression-complexity-based distance measures to analyze genomic sequences, generate phylogenetic trees, and train machine learning models using complexity values as features.

#Authors
Amanchi Prem Kumar
Munagala Naga Venkata Trinath Sai
Athira Panicker
Karthi Balasubramanian
Affiliation: Department of Electronics and Communication Engineering, Amrita School of Engineering, Coimbatore, Amrita Vishwa Vidyapeetham, India.

Nithin Nagaraj
Affiliation: Consciousness Studies Programme, National Institute of Advanced Studies, Bengaluru, Karnataka, India.

#Setup Instructions
These instructions will guide you in setting up the GENOME_ANALYSES project on your local machine.

Place the downloaded folder GENOME_ANALYSES in your MATLAB working directory and add it to the path (right-click on the GENOME_ANALYSES folder and select "Add to Path > Folders and Subfolders").
Change the MATLAB current working directory to ~/MATLAB/GENOME_ANALYSES/ (simply open this directory in MATLAB).
Running the Software
Classification
To run the classification program:

Open the workspace directory GENOME_ANALYSES.
Click the "Run" button in the MATLAB editor.
Alternatively, you can run it from the command window by typing Classification_main (without quotes) and pressing ENTER.

#Phylogeny Tree Generation
To generate the phylogeny tree:

Open the workspace directory GENOME_ANALYSES.
Click the "Run" button in the MATLAB editor.
Alternatively, you can run it from the command window by typing phylogeny_main (without quotes) and pressing ENTER.

#Common Instructions for Classification
Classification_main.m is the main program that you need to run for obtaining classification results.
The complete list of provided datasets used in the paper can be found in the "DataBase" section of the Readme file.
The default dataset is 'covid'. You can change the dataset by modifying line number 8 in Classification_main.m.
If you want to use your own datasets, make sure they consist of at least 2 clusters or subfolders in the DataBase folder. For .fasta files, create subfolders representing different clusters in the directory GENOME_ANALYSES/DataBase/TestDataSet and place the respective sequences (.fasta files) in the subfolders. For NCBI accession numbers, create a .txt file with comma-separated NCBI accession numbers for each cluster and place all .txt files (one per cluster) in the directory GENOME_ANALYSES/DataBase/TestDataSet.
You can choose the feature to be used for classification from the following options:
LZ (Lempel-Ziv complexity)
ETC (Effort-to-compress complexity)
ETC+LZ (Effort-to-compress complexity and Lempel-Ziv complexity)
Common Instructions for Phylogenetic Tree
phylogeny_main.m is the main program that you need to run for generating the phylogenetic tree.
The complete list of provided datasets used in the paper can be found in the "DataBase" section of the Readme file.
The default dataset is 'Mammals'. You can change the dataset by modifying line number 9 in phylogeny_main.m.
