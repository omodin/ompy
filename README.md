# ompy
Python module with functions for analysing results from 16S amplicon sequencing.

The following things can be done with ompy:

-Subsetting and manipulating a frequency table in several ways

-Visualizing the relative abundance of the most abundant taxa in a heatmap

-Calculating naive- and phylogenetic diversity indices based on Hill numbers. Both alpha and beta diversity can be calculated.

-Visualizing alpha diversity of samples for the full range of diversity orders (0-2)

-Visualizing dissimilarities between samples in a PCoA ordination


An example of how to use ompy follows below.

Step 1. Make sure you have python on your computer; ompy was written in python3. In windows, python distributions such as winpython or anaconda could be used. These distributions come with the python packages used in ompy. The one exception is the package python-Levenshtein, which you might have to install. This package is used in ompy to calculate phylogenetic diversity indices based on pairwise distances between sequences.

Step 2. Download the ompy_functions.py file and place it in a folder on your computer.

Step 3. In the same folder as ompy_functions.py, place the required input files. These are a frequency table, a fasta file, and a meta data file. Examples of how these files are formatted are provided, see e.g. 'example_frequency_table.txt'.
In the frequency table, the first column must be the OTU or SV names (OTU1, OTU2, etc), then follows columns representing the samples with their corresponding read counts. The final columns are the taxonomic names, starting with Kingdom (or Domain).
The fasta file contains the sequences associated with the OTUs (or SVs). The meta data file contains meta data associated with each sample. The first column should be the sample names.

Step 4. Open your python IDE (for example IDLEX in winpython), create a 'some_name.py' file and save it in the same folder as ompy.py. An example of this file is provided.
In the some_name.py file, first import the functions from ompy_functions using:

from ompy_functions import *

Then load your data as an ompy object using the function loadFiles():

obj = loadFiles(tab='example_frequency_table.txt', fasta='example_fasta.fa', meta='example_metadata.txt', sep='\t')

Now obj is a python dictionary containing five pandas dataframes.
obj['tab'] holds the frequency table with samples as columns and OTUs/SVs as rows, obj['ra'] holds the relative abundances (in %), obj['seq'] holds the sequences, obj['tax'] holds the taxonomy of each OTU/SV, and obj['meta'] holds the meta data with sample names as row indexes.

To visualize the relative abundances of the most abundant taxa in a heatmap use the function plotHeatmap():

plotHeatmap(obj, levels=['Class', 'Genus'], savename='Heatmap')

The levels input specifies the taxonomic levels to show on the y-axis in the heatmap. The OTUs/SVs are here grouped based on genus and the class level is also shown.
Specifying savename saves the figure as a .png file.

To analyze the alpha diversity, use the function plotDivAlpha():

plotDivAlpha(obj, savename='Alpha_div')

The command above will plot a figure shows naive alpha diversity (Hill numbers) as a function of the diversity order (q). savename saves the figure as a .png file.
Use several different input parameters to the function it is possible to specify which specific samples to include in the plot and how to color code the lines.
If a distance matrix is provided to the input parameter distmat, the plot will show phylogenetic (or functional) diversity numbers (see details in ompy_functions.py).

To calculate dissimilarity indices between samples use the functions naiveDivBeta() or phylDivBeta():

dis = naiveDivBeta(obj['tab'], q=1, rarefy='min')

This saves a matrix containing all naive dissimilarities of diversity order 1 in the variable dis. rarefy='min' specifies that the frequency table should be rarefied to the read count in the smallest sample.

To visualizy dissimilarities between samples in a principal coordinate analysis, use the function plotPCoA():

plotPCoA(dist=dis, meta=obj['meta'], var1='reactor', savename='PCoA')

This will plot and save a PCoA figure. The var1='reactor' input will use the column 'reactor' in the metadata to color code the data points.
