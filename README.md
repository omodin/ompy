# ompy
Python module with functions for analysing results from 16S amplicon sequencing.

The following things can be done with ompy:

-Subsetting and manipulating a frequency table in several ways

-Visualizing the relative abundance of the most abundant taxa in a heatmap

-Calculating naive- and phylogenetic diversity indices based on Hill numbers. Both alpha and beta diversity can be calculated.

-Visualizing dissimilarities between samples in a PCoA ordination


An example of how to use ompy follows below.

Step 1. Make sure you have python on your computer; ompy was written in python3. In windows, python distributions such as winpython or anaconda could be used. These distributions come with the python packages used in ompy. The one exception is the package python-Levenshtein, which you might have to install. This package is used in ompy to calculate phylogenetic diversity indices based on pairwise distances between sequences.

Step 2. Download the ompy.py file and place it in a folder on your computer.

Step 3. In the same folder as ompy.py, place the required input files. These are a frequency table, a fasta file, and a meta data file. Examples of how these files are formatted are provided, see e.g. 'example_frequency_table.txt'.
In the frequency table, the first column must be the OTU or SV names (OTU1, OTU2, etc), then follows columns representing the samples with their corresponding read counts. The final columns are the taxonomic names, starting with Kingdom (or Domain).
The fasta file contains the sequences associated with the OTUs (or SVs). The meta data file contains meta data associated with each sample. The first column should be the sample names.

Step 4. Open your python IDE (for example IDLEX in winpython), create a 'some_name.py' file and save it in the same folder as ompy.py. An example of this file is provided.

