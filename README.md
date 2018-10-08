# ompy
Python module with functions for analysing results from 16S amplicon sequencing.

The following things can be done with ompy: 
-Subsetting and manipulating a frequency table in several ways
-Visualizing the relative abundance of the most abundant taxa in a heatmap
-Calculating naive- and phylogenetic diversity indices based on Hill numbers. Both alpha and beta diversity can be calculated.
-Visualizing dissimilarities between samples in a PCoA ordination

An example of how to use ompy follows below.

Step 1. Make sure you have python on your computer. ompy was written in python3. In windows, for example, winpython could be used.
Step 2. Download the ompy.py file and place it in a folder on your computer.
Step 3. In the same folder as ompy.py, place the required input files. These are a frequency table, a fasta file, and a meta data file.
The first row in the frequency table contains the column heading. The first column is the name of the OTUs or SVs, e.g. OTU1, OTU2, etc. 

