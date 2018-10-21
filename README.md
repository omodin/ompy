# ompy
Python module with functions for analysing results from 16S amplicon sequencing.

The following things can be done:

-Subsetting and manipulating a frequency table in several ways

-Visualizing the relative abundance of the most abundant taxa in a heatmap

-Calculating diversity indices based on Hill numbers. Both alpha and beta diversity can be calculated.
The naive (or taxonomic) diversity indices do not take sequence similarity into account whereas the phylogenetic ones do.

-Visualizing alpha diversity of samples for the full range of diversity orders (0-2)

-Visualizing dissimilarities between samples in a PCoA ordination


The file ompy_functions.py contain all the functions. The file ompy_gui.py provides the source code for a graphical user interface
for some of the functions. A distribution of the program, which can be run on windows can be downloaded at omvatten.se/software.html.


The input data should be in the format specified in the example files:

example_frequency_table.txt has the OTUs or SVs as rows and the samples as columns. The first column holds the OTU name. The last columns hold the
taxonomic information if available, starting from 'Kingdom', then 'Phylum' etc.

example_fasta.fa is a fasta file with the sequences representing each OTU or SV.

example_metadata.txt holds the metadata


A few of the functions are described below. More detailed explanations are provided as annotations next to each function in the file ompy_functions.py.

The data is loaded using the function loadFiles():

obj = loadFiles(path='', tab='example_frequency_table.txt', fasta='example_fasta.fa', meta='example_metadata.txt', sep='\t')

Now obj is a python dictionary containing five pandas dataframes.
obj['tab'] holds the frequency table with samples as columns and OTUs/SVs as rows, obj['ra'] holds the relative abundances (in %), obj['seq'] holds the sequences, obj['tax'] holds the taxonomy of each OTU/SV, and obj['meta'] holds the meta data with sample names as row indexes.

To visualize the relative abundances of the most abundant taxa in a heatmap use the function plotHeatmap():

plotHeatmap(obj, levels=['Class', 'Genus'], savename='Heatmap')

The levels input specifies the taxonomic levels to show on the y-axis in the heatmap. The OTUs/SVs are here grouped based on genus and the class level is also shown.
Specifying savename saves the figure as a .png file.

To analyze the alpha diversity, use the function plotDivAlpha():

plotDivAlpha(obj, savename='Alpha_div')

The command above will plot a figure shows taxonomic alpha diversity (Hill numbers) as a function of the diversity order (q). savename saves the figure as a .png file.
Use several different input parameters to the function it is possible to specify which specific samples to include in the plot and how to color code the lines.
If a distance matrix is provided to the input parameter distmat, the plot will show phylogenetic diversity numbers.

To calculate dissimilarity indices between samples use the functions naiveDivBeta() or phylDivBeta():

dis = naiveDivBeta(obj['tab'], q=1, rarefy='min')

This saves a matrix containing all taxonomic dissimilarities of diversity order 1 in the variable dis. rarefy='min' specifies that the frequency table should be rarefied to the read count in the smallest sample.

To visualizy dissimilarities between samples in a principal coordinate analysis, use the function plotPCoA():

plotPCoA(dist=dis, meta=obj['meta'], var1='reactor', savename='PCoA')

This will plot and save a PCoA figure. The var1='reactor' input will use the column 'reactor' in the metadata to color code the data points.
