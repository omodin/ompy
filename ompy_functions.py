import pandas as pd
import numpy as np
import Levenshtein as Lv
import math
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

pd.options.mode.chained_assignment = None  # default='warn'

# Returns color- or marker lists to use in figures
# type is 'colors' or 'markers'. If plot=True a figure with available options will be shown else a list is returned
def get_colors_markers(type='colors', plot=False):
    # Sort colors by hue, saturation, value and name.
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name) for name, color in colors.items())
    color_names = [name for hsv, name in by_hsv]

    if type == 'colors' and not plot:
        numberlist = [128, 24, 38, 79, 146, 49, 152, 117, 58, 80, 119, 20]
        outputlist = []
        for i in numberlist:
            outputlist.append(color_names[i])
        return outputlist

    elif type == 'colors' and plot:
        n = len(color_names)
        ncols = 4
        nrows = n // ncols

        fig, ax = plt.subplots(figsize=(12, 10))

        # Get height and width
        X, Y = fig.get_dpi() * fig.get_size_inches()
        h = Y / (nrows + 1)
        w = X / ncols

        for i, name in enumerate(color_names):
            row = i % nrows
            col = i // nrows
            y = Y - (row * h) - h

            xi_line = w * (col + 0.05)
            xf_line = w * (col + 0.25)
            xi_text = w * (col + 0.3)

            ax.text(xi_text, y, str(i)+':'+name, fontsize=(h * 0.6),
                    horizontalalignment='left',
                    verticalalignment='center')
            ax.hlines(y + h * 0.1, xi_line, xf_line,
                      color=colors[name], linewidth=(h * 0.8))

        ax.set_xlim(0, X)
        ax.set_ylim(0, Y)
        ax.set_axis_off()

        fig.subplots_adjust(left=0, right=1,
                            top=1, bottom=0,
                            hspace=0, wspace=0)
        plt.show()

    elif type == 'markers' and not plot:
        return ['o', 's', 'v', 'X', '.', '*', 'P', 'D', '<', ',', '^', '>', '1', '2', '3', '4', '8', 'h', 'H', '+']

    elif type == 'markers' and plot:
        mlist = ['o', 's', 'v', 'X', '*', 'P', 'D', '<', '1', '^', '2', '>', '3', '4', '.']
        for i in range(len(mlist)):
            plt.scatter(i+1, i+1, marker=mlist[i], s=30)
            plt.text(i+1, i+1.5, mlist[i])
        plt.show()

# ------------------------------------------
# FUNCTIONS FOR LOADING AND SAVING DATA FILES

# Returns dictionary object containing data as pandas dataframes ('tab', 'ra', 'tax', 'seq', and 'meta')
# path is path to input files
# tab is frequency table is frequency table. SV names should be in first column, taxa should be in the final columns and start with Kingdom or Domain
# fasta is fasta file with sequences. SV names should correspond to those in tab
# meta is meta data
# sep specifies separator used in input files e.g. ',' or '\t'
def loadFiles(path='', tab='None', fasta='None', meta='None', sep=','):  # Import file and convert them to suitable format
    #Prepare tab and tax
    if tab != 'None':
        # Read count table with taxa information
        readtab = pd.read_csv(path + tab, sep=sep, header=0, index_col=0)

        #Check if taxa information is in table
        taxaavailable = 1
        if 'Kingdom' in readtab.columns:
            taxpos = readtab.columns.get_loc('Kingdom')
        elif 'Domain' in readtab.columns:
            taxpos = readtab.columns.get_loc('Domain')
        else:
            taxpos = len(readtab.columns)
            taxaavailable = 0

        readtab.index.name = 'SV'
        ctab = readtab.iloc[:, :taxpos]
        ratab = 100 * ctab / ctab.sum()

        ctab = ctab.sort_index()
        ratab = ratab.sort_index()

        # Prepare taxa dataframe
        if taxaavailable == 1:
            taxtab = readtab.iloc[:, taxpos:]
            taxtab = taxtab.sort_index()
            for name in taxtab.columns.tolist()[:-1]: #Remove items containing unknowns
                taxtab[name][taxtab[name].str.contains('nknown', na=False)] = np.nan
                taxtab[name][taxtab[name].str.contains('uncult', na=False)] = np.nan

    # Read fasta file with SV sequences
    if fasta != 'None':
        with open(path + fasta, 'r') as f:
            readfasta = f.readlines()
        fastalist = []
        seq = ''
        name = ''
        for i in range(len(readfasta)):
            if readfasta[i][0] == '>':
                name = readfasta[i][1:].strip()
                seq = ''
            elif i == len(readfasta) - 1 or readfasta[i + 1][0] == '>':
                seq = seq + readfasta[i].strip()
                fastalist.append([name, seq])
            else:
                seq = seq + readfasta[i].strip()
        # Correct fasta list based on SVs actually in count table (some might not be represented)
        if tab != 'None':
            tabSVs = list(ctab.index)
            corrfastalist = []
            for i in range(len(fastalist)):
                if fastalist[i][0] in tabSVs:
                    corrfastalist.append(fastalist[i])
            seqtab = pd.DataFrame(corrfastalist, columns=['SV', 'seq'])
        else:
            seqtab = pd.DataFrame(fastalist, columns=['SV', 'seq'])

        seqtab = seqtab.set_index('SV')
        seqtab = seqtab.sort_index()

    # Read meta data
    if meta != 'None':
        readmeta = pd.read_csv(path + meta, sep=sep, header=0, index_col=0)
    # Go through metadata and remove lines not in tab
    if meta != 'None' and tab != 'None':
        for ix in readmeta.index:
            if ix not in ctab.columns:
                readmeta = readmeta.drop([ix])

    # Return dictionary object with all dataframes
    out = {}
    if tab != 'None':
        out['tab'] = ctab
        out['ra'] = ratab
    if tab != 'None' and taxaavailable == 1:
        out['tax'] = taxtab
    if fasta != 'None':
        out['seq'] = seqtab
    if meta != 'None':
        out['meta'] = readmeta
    return out

# Outputs frequency table, fasta file and meta data from an object
# obj is object to be returned, path is path to folder where files are to be saved
def returnFiles(obj, path, sep=','):  # Saves files in the same format as they were loaded
    # Return taxa-count table
    if 'tab' in obj and 'tax' in obj:
        tab = obj['tab']; tax = obj['tax']
        tab_tax = pd.concat([tab, tax], axis=1)
        tab_tax.to_csv(path + 'output_CountTaxaTable.csv', sep=sep)
    elif 'tab' in obj:
        tab = obj['tab']
        tab.to_csv(path+'output_CountTable.csv', sep=sep)
    else:
        print('No tab and tax')

    if 'meta' in obj:
        meta = obj['meta']
        meta.to_csv(path + 'output_Meta.csv', sep=sep)
    else:
        print('No meta')

    if 'seq' in obj:
        seq = obj['seq']
        fasta = []
        for s in seq.index:
            fasta.append('>' + s + '\n')
            fasta.append(seq.loc[s, 'seq'] + '\n')
        with open(path + 'output_SVs.fa', 'w') as f:
            for i in fasta:
                f.write(i)
    else:
        print('No seq')
    print('Files saved')

# Returns some information about an object, e.g. number of samples, reads, headings in meta data etc.
def getStats(obj):
    tab = obj['tab']
    print('Total number of samples= ', len(tab.columns))
    print('Total number of SVs= ', len(tab.index))
    print('Total counts= ', sum(tab.sum()))
    print('Minimum number of counts in a sample= ', min(tab.sum()))
    print('Titles and first line in meta data')
    print(list(obj['meta'].columns))
    print(list(obj['meta'].iloc[1, :]))


# ----------------------------
# FUNCTIONS FOR SUBSETTING DATA

# Subsets samples based on metadata
# var is the column heading in metadata used to subset samples
# slist is a list of names in meta data column which specify samples to keep
# if keep0 is false, all SVs with 0 counts after the subsetting will be discarded from the data
def subsetSamples(obj, var, slist, keep0=False):
    if 'meta' not in obj.keys():
        print('Metadata missing')
        return 0

    meta = obj['meta']
    meta = meta[meta[var].isin(slist)]
    tab = obj['tab']
    tab = tab[meta.index]
    ra = obj['ra']
    ra = ra[meta.index]
    if 'seq' in obj.keys():
        seq = obj['seq']
    if 'tax' in obj.keys():
        tax = obj['tax']
    out = {}  # Dictionary object to return, dataframes and dictionaries

    # Remove SV with zero count
    if not keep0:
        tab['sum'] = tab.sum(axis=1)
        tab2 = tab[tab['sum'] > 0]
        keepSVs = list(tab2.index)
        tab2 = tab2.drop('sum', axis=1)
        out['tab'] = tab2
        ra2 = ra.loc[keepSVs, :]
        out['ra'] = ra2
        if 'seq' in obj.keys():
            seq2 = seq.loc[keepSVs, :]
            out['seq'] = seq2
        if 'tax' in obj.keys():
            tax2 = tax.loc[keepSVs, :]
            out['tax'] = tax2
        out['meta'] = meta
    else:
        out['meta'] = meta
        out['tab'] = tab
        out['ra'] = ra
        if 'seq' in obj.keys():
            out['seq'] = seq
        if 'tax' in obj.keys():
            out['tax'] = tax
    return out

# Subsets object based on list of SVs to keep
def subsetSVs(obj, svlist):
    out = {}
    tab = obj['tab']
    tab = tab.loc[svlist, :]
    out['tab'] = tab
    if 'ra' in obj:
        ra = obj['ra']
        ra = ra.loc[svlist, :]
        out['ra'] = ra
    if 'tax' in obj:
        tax = obj['tax']
        tax = tax.loc[svlist, :]
        out['tax'] = tax
    if 'seq' in obj:
        seq = obj['seq']
        seq = seq.loc[svlist, :]
        out['seq'] = seq
    if 'meta' in obj:
        out['meta'] = obj['meta']
    return out

# Subsets object to the most abundant SVs calculated based on the sum of the relative abundances in each sample
# number specifies the number of SVs to keep
def subsetTopSVs(obj, number=25):
    out = {}
    tab = obj['tab']
    ra = tab/tab.sum()
    ra['sum'] = ra.sum(axis=1)
    ra = ra.sort_values(by='sum', ascending=False)
    svlist = ra.index[:number]
    tab2 = tab.loc[svlist, :]
    out['tab'] = tab2
    if 'ra' in obj:
        ra2 = obj['ra']
        ra2 = ra2.loc[svlist, :]
        out['ra'] = ra
    if 'tax' in obj:
        tax = obj['tax']
        tax = tax.loc[svlist, :]
        out['tax'] = tax
    if 'seq' in obj:
        seq = obj['seq']
        seq = seq.loc[svlist, :]
        out['seq'] = seq
    if 'meta' in obj:
        out['meta'] = obj['meta']
    return out

# Subset object based on text in taxonomic names
# subsetLevels is list taxonomic levels searched for text patterns, e.g. ['Family', 'Genus']
# subsetPatterns is list of text to search for, e.g. ['Nitrosom', 'Brochadia']
def subsetTextPatterns(obj, subsetLevels=[], subsetPatterns=[]):
    if len(subsetLevels) == 0 or len(subsetPatterns) == 0:
        print('No taxlevels or pattern')
        return 0
    else:
        taxPatternCheck = obj['tax'].applymap(str)
        keepIX = []
        for col in subsetLevels:
            for ix in taxPatternCheck.index:
                for ptrn in subsetPatterns:
                    if ptrn in taxPatternCheck.loc[ix, col] and ix not in keepIX:
                        keepIX.append(ix)
    out = {}
    if 'tab' in obj.keys():
        out['tab'] = obj['tab'].loc[keepIX]
    if 'ra' in obj.keys():
        out['ra'] = obj['ra'].loc[keepIX]
    if 'tax' in obj.keys():
        out['tax'] = obj['tax'].loc[keepIX]
    if 'seq' in obj.keys():
        out['seq'] = obj['seq'].loc[keepIX]
    if 'meta' in obj.keys():
        out['meta'] = obj['meta']
    return out

# Merges samples based on information in meta data
# var is the column heading in metadata used to merge samples. The counts for all samples with the same text in var column will be merged.
# slist is a list of names in meta data column which specify samples to keep. If slist='None' (default), the whole meta data column is used
# if keep0 is false, all SVs with 0 counts after the subsetting will be discarded from the data
def mergeSamples(obj, var='None', slist='None', keep0=False):
    if 'meta' not in obj.keys():
        print('Metadata missing')
        return 0
    if var != 'None' and slist == 'None':
        slist = obj['meta'][var]

    tabdi = {}
    radi = {}  # Temp dict that holds sum for each type in slist
    for smp in slist:
        tempobj = subsetSamples(obj, var, [smp], keep0=True)
        tab = tempobj['tab']
        tab['sum'] = tab.sum(axis=1)
        tab['ra'] = 100 * tab['sum'] / tab['sum'].sum()
        tabdi[smp] = tab.loc[:, 'sum']
        radi[smp] = tab.loc[:, 'ra']
    temptab = pd.DataFrame(tabdi, index=obj['tab'].index)
    tempra = pd.DataFrame(radi, index=obj['tab'].index)

    out = {}
    if keep0 == False:  ## Remove SV with zero count
        temptab['sum'] = temptab.sum(axis=1)
        tab2 = temptab[temptab['sum'] > 0]
        keepSVs = list(tab2.index)
        tab2 = tab2.drop('sum', axis=1)
        ra2 = tempra.loc[keepSVs, :]
        if 'seq' in obj.keys():
            seq = obj['seq']
            seq2 = seq.loc[keepSVs, :]
            out['seq'] = seq2
        if 'tax' in obj.keys():
            tax = obj['tax']
            tax2 = tax.loc[keepSVs, :]
            out['tax'] = tax2
    else:
        tab2 = temptab
        ra2 = tempra
        if 'seq' in obj.keys():
            out['seq'] = obj['seq']
        if 'tax' in obj.keys():
            out['tax'] = obj['tax']
    out['tab'] = tab2
    out['ra'] = ra2

    meta = obj['meta']
    meta2 = meta.groupby(var).first()
    meta2[var] = meta2.index
    out['meta'] = meta2
    return out

# Rarefies frequency table to a specific number of reads per sample
# if depth = 'min', the minimum number of reads in a sample is used
# seed sets a random state for reproducible results
# The function is  similar to rarefaction without replacement
def rarefy1(tab, depth='min', seed='None'):
    tab = tab.applymap(int) #Make sure table elements are integers
    if depth == 'min':  # Define read depth
        reads = min(tab.sum())
    else:
        reads = depth

    samples = tab.columns.tolist()
    rtab = tab.copy()
    for smp in samples:
        smpsum = sum(tab[smp])
        if smpsum < reads:  # Remove sample if sum of reads less than read depth
            rtab = rtab.drop(smp, axis=1)
            continue

        frac = tab.loc[:, smp] * reads / smpsum #Calculate rel abund
        avrundad = frac.apply(math.floor) #Round down
        addNR = int(reads - sum(avrundad)) #This many reads must be added

        if addNR >= 1:
            diffs = frac - avrundad
            if seed == 'None':
                addSV = tab[smp].sample(n=addNR, replace=False, weights=diffs).index.tolist()
            else:
                addSV = tab[smp].sample(n=addNR, replace=False, weights=diffs, random_state=seed).index.tolist()
            avrundad[addSV] = avrundad[addSV] + 1
        rtab[smp] = avrundad
    return rtab

# Rarefies frequency table to a specific number of reads per sample
# if depth = 'min', the minimum number of reads in a sample is used
# seed sets a random state for reproducible results
# This is rarefaction with replacement
def rarefy2(tab, depth='min', seed='None'):
    tab = tab.fillna(0)
    if seed != 'None':
        prng = np.random.RandomState(seed) # reproducible results
    noccur = tab.sum()
    nvar = len(tab.index) # number of SVs

    ## Set read depth
    if depth == 'min':
        depth = int(np.min(noccur))
    else:
        depth = int(depth)

    rftab = tab.copy()
    for s in tab.columns: # for each sample
        if tab[s].sum() < depth:
            rftab = rftab.drop(s, axis=1)
            continue
        else:
            p = tab[s]/tab[s].sum()
            if seed != 'None':
                choice = prng.choice(nvar, depth, p=p)
            else:
                choice = np.random.choice(nvar, depth, p=p)
            rftab[s] = np.bincount(choice, minlength=nvar)
    return rftab

# ------------------------------------------
# FUNCTIONS FOR VISUALISING TAXA IN HEATMAP

# Groups SVs based on taxa, returns object with grouped sequences
# levels specifies taxonomic levels to use in the grouping
# nameType specifies the abbreviation to be used for unclassified sequences
def groupbyTaxa(obj, levels=['Phylum', 'Genus'], nameType='SV'):
    # Clean up tax data frame
    tax = obj['tax']
    tax = tax.fillna(0)
    taxSV = tax.copy() #df to hold nameTypes in undefined
    taxNames = tax.copy() #df to hold lowest known taxaname in undefined

    # Check which OTU or SV name is used in the index
    indexlist = tax.index.tolist()
    if indexlist[0][:3] in ['Otu', 'OTU', 'ASV', 'ESV']:
        currentname = indexlist[0][:3]
        startpos = 3
    elif indexlist[0][:2] in ['SV']:
        currentname = indexlist[0][:2]
        startpos = 2
    else:
        print('Error in groupbyTaxa, SV/OTU name not known')
        return 0

    # If incorrect name is in tax, make column with correct name
    if nameType != currentname:
        newnames = []
        for i in range(len(indexlist)):
            newnames.append(nameType+indexlist[i][startpos:])
        indexlist = newnames

    # Put the SV/OTU name in all empty spots in taxSV
    for col in range(len(taxSV.columns)):
        for row in range(len(taxSV.index)):
            if taxSV.iloc[row, col] == 0:
                taxSV.iloc[row, col] = indexlist[row]
    taxSV[nameType] = indexlist

    # Change all 0 in tax to lowest determined taxa level in taxNames
    taxanameslist = taxNames.columns.tolist() #List with Kingdom, Phylum .. SV
    for s_nr in range(1, len(taxanameslist)):
        s0 = taxanameslist[s_nr-1]
        s1 = taxanameslist[s_nr]
        taxNames[s1][tax[s1] == 0] = taxNames[s0][tax[s1] == 0]

    # Create names to use in output
    if len(levels) == 1:
        tax['Name'] = taxSV[levels[0]]
        for ix in tax.index:
            if tax.loc[ix, levels[0]] == 0:
                tax.loc[ix, 'Name'] = taxNames.loc[ix, levels[0]] + '; ' + tax.loc[ix, 'Name']
    elif len(levels) == 2:
        tax['Name'] = taxNames[levels[0]]+'; '+taxSV[levels[1]]
    else:
        print('Error in GroupbyTaxa, levels should be a list with 1 or 2 items')
        return 0

    #Grouby Name and return object
    out = {}
    if 'tab' in obj.keys():
        tab = obj['tab']
        tab['Name'] = tax['Name']
        tab = tab.set_index(['Name'])
        tab = tab.groupby(tab.index).sum()
        out['tab'] = tab
    if 'ra' in obj.keys():
        ra = obj['ra']
        ra['Name'] = tax['Name']
        ra = ra.set_index('Name')
        ra = ra.groupby(ra.index).sum()
        out['ra'] = ra
    if 'meta' in obj.keys():
        out['meta'] = obj['meta']
    return out

# Plots heatmap
# var specifies heading in meta data used to merge samples
# levels specifies taxonomic levels used in y axis
# subsetLevels and subsetPatterns refer to subsetTextPatters function which can be used to filter results
# order refers to heading in meta data used to order samples
# numberToPlot refers to the number of taxa with highest abundance to include in the heatmap
# method refers to the method used to define the taxa with highest abundance, 'max_sample' is max relative abundance in a sample,
# 'mean_all' is the mean relative abundance across all samples
# nameType is nameType in groupbyTaxa function
# if labels=True, include relative abundance values in heatmap, if False they are not included
# labelSize is the size of the relative abundance lables in the heatmap
# fontSize is the size of the axis text
# savename is the name (also include path) of the saved png file, if 'None' no figure is saved
def plotHeatmap(obj, var='None', levels=['Phylum', 'Genus'], subsetLevels='None', subsetPatterns='None',
                order='None', numberToPlot=20, method='max_sample', nameType='SV', labels=True, labelSize=1, fontSize=15, savename='None'):
    #Merge samples based on var
    if var != 'None':
        merged_obj = mergeSamples(obj, var=var)
    else:
        merged_obj = obj.copy()

    #Calculate relative abundances and store in df ra
    tab = merged_obj['tab']
    ra = 100*tab/tab.sum()
    merged_obj['ra'] = ra

    ## Make sure samples are in the right order in meta data
    if order != 'None':
        md = merged_obj['meta']
        md[order] = md[order].astype(float)
        md = md.sort_values(by=order)
        logiclist = []
        if var != 'None':
            [logiclist.append(item) for item in md[var] if item not in logiclist]
        else:
            [logiclist.append(item) for item in md.index if item not in logiclist]
        merged_obj['meta'] = md

    ## Subset based on pattern
    if subsetLevels != 'None' and isinstance(subsetLevels, list) and isinstance(subsetPatterns, list):
        subset_obj = subsetTextPatterns(merged_obj, subsetLevels, subsetPatterns)
        merged_obj = subset_obj

    ## Groupby taxa
    taxa_obj = groupbyTaxa(merged_obj, levels=levels, nameType=nameType)
    ra = taxa_obj['ra']; table = ra.copy()

    # Subset for top taxa
    if method == 'max_sample':
        ra['max'] = ra.max(axis=1)
        ra = ra.sort_values(by=['max'], ascending=False)
        retain = ra.index.tolist()[:numberToPlot]

    elif method == 'mean_all':
        ra['mean'] = ra.mean(axis=1)
        ra = ra.sort_values(by=['mean'], ascending=False)
        retain = ra.index.tolist()[:numberToPlot]

    table = table.loc[retain]

    if order != 'None':
        table = table.loc[:, logiclist]

    # Change to italics in table labels
    taxa_list = table.index.tolist()
    new_taxa_list = []
    for n in taxa_list:
        if ';' in n: #Check if there are two taxa names
            splitname = n.split(';')
            splitname1 = splitname[0].split('__')
            newname1 = splitname1[0]+'__'+'$\it{'+splitname1[1]+'}$'
            if '__' in splitname[1]:
                splitname2 = splitname[1].split('__')
                newname2 = splitname2[0]+'__'+'$\it{'+splitname2[1]+'}$'
            else:
                newname2 = splitname[1]
            newname = newname1+';'+newname2
        else: #If there is only one taxa name
            if '__' in n:
                splitname = n.split('__')
                newname = splitname[0]+'__'+'$\it{'+splitname[1]+'}$'
            else:
                newname = n
        new_taxa_list.append(newname)
    table = pd.DataFrame(table.values, index=new_taxa_list, columns=table.columns)

    # Print heatmap
    table['avg'] = table.mean(axis=1)
    table = table.sort_values(by=['avg'], ascending=True)
    table = table.drop(['avg'], axis=1)

    #Fix datalabels
    if labels:
        labelvalues = table.copy()
        for r in table.index:
            for c in table.columns:
                value = float(table.loc[r, c])
                if value < 0.1 and value > 0:
                    labelvalues.loc[r, c] = '<0.1'
                elif value < 10 and value >= 0.1:
                    labelvalues.loc[r, c] = round(value, 1)
                elif value > 99:
                    labelvalues.loc[r, c] = int(99)
                elif value >= 10:
                    labelvalues.loc[r, c] = round(value)
                else:
                    labelvalues.loc[r, c] = 0
        labelvalues = labelvalues.applymap(str)

    #Plot
    fig, ax = plt.subplots(figsize=(14, 10))
    sns.set(font_scale=labelSize)
    if labels:
        sns.heatmap(table, annot=labelvalues, fmt='', cmap='Reds', linewidths=0.5, robust=True, cbar=False, ax=ax)
    else:
        sns.heatmap(table, cmap='Reds', linewidths=0.5, robust=True, cbar=False, ax=ax)
    plt.xticks(rotation=90)
    plt.setp(ax.get_xticklabels(), fontsize=fontSize)
    ax.set_ylabel('')
    plt.yticks(rotation=0)
    plt.setp(ax.get_yticklabels(), fontsize=fontSize)
    plt.tight_layout()
    if savename != 'None':
        plt.savefig(savename)
    plt.show()

# ------------------------------------------
# FUNCTIONS FOR ANALYSING DIVERSITY

# Returns matrix for pairwise distances between SVs based on Levenshtein/Hamming distances (uses Levenshtein package)
# Saves results as pickle file at location specified in savename
# Output is needed as input for phylogenetic diversity index calculations
def phylDistMat(seq, pickleOrCsv='pickle', savename='PhylDistMat'):
    svnames = list(seq.index)
    df = pd.DataFrame(0, index=svnames, columns=svnames)
    for i in range(len(svnames) - 1):
        if i % 10 == 0:
            print(i, ' out of ', len(svnames))
        for j in range(i + 1, len(svnames)):
            n1 = svnames[i]
            n2 = svnames[j]
            s1 = seq.loc[n1, 'seq']
            s2 = seq.loc[n2, 'seq']

            if len(s1) == len(s2):
                dist = Lv.hamming(s1, s2) / len(s1)
            else:
                maxlen = max(len(s1), len(s2))
                dist = Lv.distance(s1, s2) / maxlen

            df.loc[n1, n2] = dist; df.loc[n2, n1] = dist
    if pickleOrCsv == 'pickle':
        with open(savename+'.pickle', 'wb') as handle:
            pickle.dump(df, handle)
    else:
        df.to_csv(savename+'.csv')
    print('Finished printing Phyldistmat as pickle or csv file')

# Returns Rao's quadratic entropy, sum(sum(dij*pi*pj))
# Function used in Chiu's phylogenetic diversity functions
def raoQ(tab, distmat):
    ra = tab / tab.sum()
    outdf = pd.Series(0, index=ra.columns)
    svlist = list(ra.index)
    distmat = distmat.loc[svlist, svlist]
    for smp in ra.columns:
        ra2mat = pd.DataFrame(np.outer(ra.loc[:, smp].values, ra.loc[:, smp].values), index=ra.index, columns=ra.index)
        rao_mat = ra2mat.mul(distmat)
        Qvalue = sum(rao_mat.sum())
        outdf.loc[smp] = Qvalue
    return outdf

# Converts beta value to distances, specify q and type associated with the beta (assuming pairwise)
# Used in beta dissimilarity calculations
def beta2Dist(beta, q=1, divType='naive'):
    beta = beta.applymap(float)
    dist = beta.copy()
    mask = beta > 0
    if divType == 'naive':
        if q == 1:
            dist[mask] = 1 - (math.log(2) - np.log(beta[mask])) / math.log(2)
        else:
            dist[mask] = 1 - ((1 / beta[mask]) ** (q - 1) - 0.5 ** (q - 1)) / (1 - 0.5 ** (q - 1))

    elif divType == 'phyl':
        if q == 1:
            dist[mask] = np.log(beta[mask]) / (2 * math.log(2))
        else:
            dist[mask] = 1 - ((2 ** (2 * (1 - q)) - beta[mask].pow(1 - q)) / (2 ** (2 * (1 - q)) - 1))
    return dist

# Returns naive alpha diversity of order q for all samples
def naiveDivAlpha(tab, q=1, rarefy='None'):
    if rarefy != 'None':
        raretab = rarefy1(tab, rarefy)
        ra = raretab / raretab.sum()
    else:
        ra = tab / tab.sum()

    if q == 0:
        mask = tab > 0
        Hillvalues = tab[mask].count()
        return Hillvalues
    elif q == 1:
        mask = (ra != 0)
        raLn = ra[mask] * np.log(ra[mask])
        Hillvalues = np.exp(-raLn.sum())
        return Hillvalues
    else:
        mask = (ra != 0)
        ra = ra[mask]
        rapow = ra[mask].pow(q)
        rapowsum = rapow.sum()
        Hillvalues = rapowsum.pow(1 / (1 - q))
        return Hillvalues

# Returns matrix of paiwise dissimilarities of order q
# Based on local overlaps as defined in Chao et al. 2014
def naiveDivBeta(tab, q=1, rarefy='None', dis=True):
    if rarefy != 'None':
        raretab = rarefy1(tab, rarefy)
        ra = raretab / raretab.sum()
    else:
        ra = tab / tab.sum()

    smplist = ra.columns
    outdf = pd.DataFrame(0, index=smplist, columns=smplist)
    for smp1nr in range(len(smplist) - 1):
        for smp2nr in range(smp1nr + 1, len(smplist)):
            smp1 = smplist[smp1nr]
            smp2 = smplist[smp2nr]

            if q == 0:
                mask1 = ra[smp1] != 0
                alpha1 = ra[smp1][mask1].count()
                mask2 = ra[smp2] != 0
                alpha2 = ra[smp2][mask2].count()
                alphavalue = (alpha1 + alpha2) / 2

                ra_sum = ra[[smp1, smp2]].sum(axis=1)
                maskg = (ra_sum != 0)
                gammavalue = ra_sum[maskg].count()

                betavalue = gammavalue / alphavalue
                outdf.loc[smp1, smp2] = betavalue
                outdf.loc[smp2, smp1] = betavalue
            elif q == 1:
                mask1 = ra[smp1] != 0
                raLn1 = ra[smp1][mask1] * np.log(ra[smp1][mask1])
                raLn1 = raLn1.sum()
                mask2 = ra[smp2] != 0
                raLn2 = ra[smp2][mask2] * np.log(ra[smp2][mask2])
                raLn2 = raLn2.sum()
                alphavalue = math.exp(-0.5 * raLn1 - 0.5 * raLn2)

                ra_mean = ra[[smp1, smp2]].mean(axis=1)
                maskg = (ra_mean != 0)
                raLng = ra_mean[maskg] * np.log(ra_mean[maskg])
                raLng = raLng.sum()
                gammavalue = math.exp(-raLng)

                betavalue = gammavalue / alphavalue
                outdf.loc[smp1, smp2] = betavalue
                outdf.loc[smp2, smp1] = betavalue
            else:
                mask1 = ra[smp1] != 0
                ra1 = ra[smp1][mask1]
                ra1pow = ra1.pow(q)
                ra1powsum = ra1pow.sum()
                mask2 = ra[smp2] != 0
                ra2 = ra[smp2][mask2]
                ra2pow = ra2.pow(q)
                ra2powsum = ra2pow.sum()
                alphavalue = (0.5 * ra1powsum + 0.5 * ra2powsum) ** (1 / (1 - q))

                ra_mean = ra[[smp1, smp2]].mean(axis=1)
                maskg = (ra_mean != 0)
                rag = ra_mean[maskg]
                ragpow = rag.pow(q)
                ragpowsum = ragpow.sum()
                gammavalue = ragpowsum ** (1 / (1 - q))

                betavalue = gammavalue / alphavalue
                outdf.loc[smp1, smp2] = betavalue
                outdf.loc[smp2, smp1] = betavalue

    if dis:
        dist = beta2Dist(beta=outdf, q=q, divType='naive')
        return dist
    else:
        return outdf

# Calculate matrix of pairwise Bray-Curtis dissimilarities
def bray(tab):
    ra = tab / tab.sum()
    smplist = list(tab.columns)
    outdf = pd.DataFrame(0, index=smplist, columns=smplist)
    for smp1nr in range(len(smplist) - 1):
        smp1 = smplist[smp1nr]
        for smp2nr in range(smp1nr + 1, len(smplist)):
            smp2 = smplist[smp2nr]
            brayvalue = 1 - (ra.loc[:, [smp1, smp2]].min(axis=1).sum())
            outdf.loc[smp1, smp2] = brayvalue;
            outdf.loc[smp2, smp1] = brayvalue
    return outdf

# Returns phylogenetic alpha diversity of order q for all samples
# FD as in Chiu et al. Plos One, 2014
def phylDivAlpha(tab, distmat, q=0, rarefy='None'):
    if rarefy != 'None':
        raretab = rarefy1(tab, rarefy)
        ra = raretab / raretab.sum()
    else:
        ra = tab / tab.sum()

    outdf = pd.Series(0, index=ra.columns)
    svlist = list(ra.index)
    distmat = distmat.loc[svlist, svlist]
    Qframe = raoQ(tab, distmat)
    if q == 0:
        for smp in tab.columns:
            ra2mat = pd.DataFrame(np.outer(ra.loc[:, smp].values, ra.loc[:, smp].values), index=ra.index,
                                  columns=ra.index)
            mask = ra2mat > 0
            ra2mat[mask] = 1
            dQmat = distmat.mul(1 / Qframe.loc[smp])
            ra2dq = (ra2mat.mul(dQmat))
            Chiuvalue = pow(sum(ra2dq.sum()), 1 / (2 * (1 - q)))
            outdf.loc[smp] = Chiuvalue
    elif q == 1:
        for smp in tab.columns:
            ra2mat = pd.DataFrame(np.outer(ra.loc[:, smp].values, ra.loc[:, smp].values), index=ra.index,
                                  columns=ra.index)
            mask = (ra2mat != 0)
            ra2Lnmat = np.log(ra2mat[mask])
            ra2ochLn = ra2mat.mul(ra2Lnmat)
            dQmat = distmat.mul(1 / Qframe.loc[smp])
            dQ_ra2_Ln = dQmat.mul(ra2ochLn)
            Chiuvalue = math.exp(-0.5 * sum(dQ_ra2_Ln.sum()))
            outdf.loc[smp] = Chiuvalue
    else:
        for smp in tab.columns:
            ra2mat = pd.DataFrame(np.outer(ra.loc[:, smp].values, ra.loc[:, smp].values), index=ra.index,
                                  columns=ra.index)
            ra2matq = ra2mat.pow(q)
            dQmat = distmat.mul(1 / Qframe.loc[smp])
            ra2dq = (ra2matq.mul(dQmat))
            Chiuvalue = pow(sum(ra2dq.sum()), 1 / (2 * (1 - q)))
            outdf.loc[smp] = Chiuvalue
    MD = outdf.mul(Qframe)
    return outdf.mul(MD)

# Returns matrix of paiwise phylogenetic dissimilarities of order q
# Based on local functional overlaps as defined in Chao et al. 2014
def phylDivBeta(tab, distmat, q=1, rarefy='None', dis=True):
    if rarefy != 'None':
        raretab = rarefy1(tab, rarefy)
        ra = raretab / raretab.sum()
    else:
        ra = tab / tab.sum()

    smplist = list(ra.columns)
    outD = pd.DataFrame(0, index=smplist, columns=smplist)

    for smp1nr in range(len(smplist) - 1):
        print('Running ', smp1nr, ' out of ', len(smplist))
        for smp2nr in range(smp1nr + 1, len(smplist)):
            smp1 = smplist[smp1nr]
            smp2 = smplist[smp2nr]

            ra12 = ra.loc[:, [smp1, smp2]]
            ra12['mean'] = ra12.mean(axis=1)
            Qvalues = raoQ(ra12, distmat)
            Qpooled = Qvalues['mean']
            dqmat = distmat.mul(1 / Qpooled)

            if q != 1:
                # Get gamma
                mask = ra12['mean'] > 0
                ra2mat = pd.DataFrame(np.outer(ra12['mean'][mask], ra12['mean'][mask]), index=ra12[mask].index,
                                      columns=ra12[mask].index)
                ra2matq = ra2mat.pow(q)
                ra2dq = ra2matq.mul(dqmat)
                Dg = pow(sum(ra2dq.sum()), 1 / (2 * (1 - q)))

                # Get alpha
                mask1 = ra12[smp1] > 0
                ra2mat = pd.DataFrame(np.outer(ra12[mask1][smp1], ra12[mask1][smp1]), index=ra12[mask1].index,
                                      columns=ra12[mask1].index)
                ra2mat = ra2mat / 4
                ra2matq = ra2mat.pow(q)
                ra2dq = ra2matq.mul(dqmat)
                asum1 = sum(ra2dq.sum())

                mask2 = ra12[smp2] > 0
                ra2mat = pd.DataFrame(np.outer(ra12[mask2][smp2], ra12[mask2][smp2]), index=ra12[mask2].index,
                                      columns=ra12[mask2].index)
                ra2mat = ra2mat / 4
                ra2matq = ra2mat.pow(q)
                ra2dq = ra2matq.mul(dqmat)
                asum2 = sum(ra2dq.sum())

                ra2mat = pd.DataFrame(np.outer(ra12[mask1][smp1], ra12[mask2][smp2]), index=ra12[mask1].index,
                                      columns=ra12[mask2].index)
                ra2mat = ra2mat / 4
                ra2matq = ra2mat.pow(q)
                ra2dq = ra2matq.mul(dqmat)
                asum12 = sum(ra2dq.sum())

                Da = 0.5 * pow((asum1 + asum2 + 2 * asum12), 1 / (2 * (1 - q)))

                # Calculate beta
                outD.loc[smp1, smp2] = Dg / Da;
                outD.loc[smp2, smp1] = Dg / Da

            else:
                # Get gamma
                mask = ra12['mean'] > 0
                ra2mat = pd.DataFrame(np.outer(ra12['mean'][mask], ra12['mean'][mask]), index=ra12[mask].index,
                                      columns=ra12[mask].index)
                ra2matq = ra2mat * np.log(ra2mat)
                ra2dq = ra2matq.mul(dqmat)
                Dg = math.exp(-0.5 * sum(ra2dq.sum()))

                # Get alpha
                mask1 = ra12[smp1] > 0
                ra2mat = pd.DataFrame(np.outer(ra12[mask1][smp1], ra12[mask1][smp1]), index=ra12[mask1].index,
                                      columns=ra12[mask1].index)
                ra2mat = ra2mat / 4
                ra2matq = ra2mat * np.log(ra2mat)
                ra2dq = ra2matq.mul(dqmat)
                asum1 = sum(ra2dq.sum())

                mask2 = ra12[smp2] > 0
                ra2mat = pd.DataFrame(np.outer(ra12[mask2][smp2], ra12[mask2][smp2]), index=ra12[mask2].index,
                                      columns=ra12[mask2].index)
                ra2mat = ra2mat / 4
                ra2matq = ra2mat * np.log(ra2mat)
                ra2dq = ra2matq.mul(dqmat)
                asum2 = sum(ra2dq.sum())

                ra2mat = pd.DataFrame(np.outer(ra12[mask1][smp1], ra12[mask2][smp2]), index=ra12[mask1].index,
                                      columns=ra12[mask2].index)
                ra2mat = ra2mat / 4
                ra2matq = ra2mat * np.log(ra2mat)
                ra2dq = ra2matq.mul(dqmat)
                asum12 = sum(ra2dq.sum())

                Da = 0.5 * math.exp(-0.5 * (asum1 + asum2 + 2 * asum12))

                # Calculate beta
                outD.loc[smp1, smp2] = Dg / Da;
                outD.loc[smp2, smp1] = Dg / Da
    outFD = outD.pow(2)

    if dis:
        return beta2Dist(beta=outFD, q=q, divType='phyl')
    else:
        return outFD

# Visualizes how alpha diversity depends on diversity order
# If distmat is specified phylogenetic alpha diversity is calculated, else naive
# rarefy specifies depth to rarefy frequency table (default is the smallest sample count, 'min')
# var refers to column heading in meta data used to color code samples
# slist is a list of samples from the var column to include (default is all)
# order refers to column heading in meta data used to order the samples
# If ylog=True, the y-axis of the plot will be logarithmic
def plotDivAlpha(obj, distmat='None', rarefy='min', var='None', slist='All', order='None', ylog=False, colorlist='None', savename='None'):
    #Pick out samples to include based on var and slist
    meta = obj['meta']
    if order != 'None':
        meta = meta.sort_values(order)

    if var == 'None':
        smplist = meta.index.tolist()
    elif slist == 'All':
        smplist = meta[var].index.tolist()
    else:
        smplist = meta.loc[slist, var].index.tolist()

    #Dataframe for holding results
    xvalues = np.arange(0, 2.01, 0.02)
    df = pd.DataFrame(0, index=xvalues, columns=smplist)

    #Put data in dataframe
    tab = obj['tab'][smplist]
    for x in xvalues:
        timeprogress = int(100*x/xvalues[-1])
        if timeprogress%10 == 0:
            print(timeprogress,'% complete')
        if distmat == 'None':
            alphadiv = naiveDivAlpha(tab, q=x, rarefy=rarefy)
        else:
            alphadiv = phylDivAlpha(tab, distmat, q=x, rarefy=rarefy)
        df.loc[x, smplist] = alphadiv

    #Plot data
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(figsize=(10, 6))

    if colorlist == 'None':
        colorlist = get_colors_markers('colors')
    else:
        colorlist = colorlist
    colorcheck = []

    for s in df.columns:
        if var != 'None':
            cv = meta.loc[s, var]
        else:
            cv = s

        if cv not in colorcheck:
            colorcheck.append(cv)
            lab = cv
        else:
            lab = '_nolegend_'

        colnr = colorcheck.index(cv)
        col = colorlist[colnr % len(colorlist)]

        if ylog:
            ax.semilogy(df.index, df[s], lw=1, color=col, label=lab)
        else:
            ax.plot(df.index, df[s], lw=1, color=col, label=lab)

    ax.set_ylabel('Diversity number ($^q$D)')
    ax.set_xlabel('Diversity order (q)')
    ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])
    ax.set_xlim(0, 2)

    plt.legend()
    plt.tight_layout()
    if savename != 'None':
        plt.savefig(savename)
    plt.show()

# Visualizes dissimilarities in PCoA plot
# dist is distance matrix and meta is meta data
# var1 is heading in meta used to color code, var2 is heading in meta used to code by marker type
# tag is heading in meta used to add labels to each point in figure
# order is heading in meta used to order samples
# title is title of figure
# colorlist specifies colorlist to use for var1; same for markerlist and var2
# savename is path and name to save png figure output
def plotPCoA(dist, meta, var1='None', var2='None', tag='None', order='None', title='', colorlist='None', markerlist='None', savename='None'):
    def get_eig(d): #Function for centering and eigen-decomposition of distance matrix
        dist2 = -0.5*(d**2)
        col_mean = dist2.mean(axis=0)
        row_mean = dist2.mean(axis=1)
        tot_mean = np.array(dist2).flatten().mean()
        dist2_cent = dist2.subtract(col_mean, axis=1)
        dist2_cent = dist2_cent.subtract(row_mean, axis=0)
        dist2_cent = dist2_cent.add(tot_mean)
        vals, vects = np.linalg.eig(dist2_cent)
        return [vals, vects]

    ev_ev = get_eig(dist)
    #Correction method for negative eigenvalues
    if min(ev_ev[0])<0: #From Legendre 1998, method 1 (derived from Lingoes 1971) chapter 9, page 502
        d2 = (dist[dist != 0]**2+2*abs(min(ev_ev[0])))**0.5
        d2 = d2.fillna(0)
        ev_ev = get_eig(d2)

    #Get proportions and coordinates
    vals = ev_ev[0].copy()
    vects = ev_ev[1]
    prop = []
    coords = []
    for i in range(2):
        maxpos = np.argmax(vals)
        prop.append(vals[maxpos]/sum(ev_ev[0]))
        coords.append((vals[maxpos]**0.5)*vects[:, maxpos])
        vals[maxpos] = 0

    # Do the plotting

    #Set axis names and make dataframe for plotting
    pc1_perc = round(100 * prop[0], 1)
    xn = 'PC1 (' + str(pc1_perc) + '%)'
    if '+' in xn:
        xn = xn[:4] + xn[-1]
    pc2_perc = round(100 * prop[1], 1)
    yn = 'PC2 (' + str(pc2_perc) + '%)'
    if '+' in yn:
        yn = yn[:4] + yn[-1]

    smplist = dist.index
    pcoadf = pd.DataFrame({xn:coords[0], yn:coords[1]}, index=smplist)

    # Combine pcoa results with meta data
    meta[xn] = pcoadf[xn]
    meta[yn] = pcoadf[yn]
    metaPlot = meta[meta[xn].notnull()]
    if order != 'None':
        meta = meta.sort_values(by=[order])

    if var1 == 'None' and var2 == 'None':
        return 'Error, no variables in input'
    if var1 != 'None': #List of names used for different colors in legend
        smpcats1 = []
        [smpcats1.append(item) for item in meta[var1] if item not in smpcats1]
    if var2 != 'None': #List of names used for different marker types in legend
        smpcats2 = []
        [smpcats2.append(item) for item in meta[var2] if item not in smpcats2]
    if tag != 'None': #List of labels placed next to points
        tagcats = []
        [tagcats.append(item) for item in meta[tag] if item not in tagcats]

    # Create figure
    if colorlist == 'None':
        colorlist = get_colors_markers('colors')
    if markerlist == 'None':
        markerlist = get_colors_markers('markers')

    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(10, 6))

    linesColor = [[], []]
    linesShape = [[], []]
    shapeTracker = []

    for i in range(len(smpcats1)):
        metaPlot_i = metaPlot[metaPlot[var1] == smpcats1[i]] #Subset metaPlot based on var1 in smpcats1

        if var2 != 'None':
            linesColor[0].append(ax.scatter([], [], label=str(smpcats1[i]), color=colorlist[i]))
            linesColor[1].append(smpcats1[i])

            jcounter = 0
            for j in range(len(smpcats2)):
                if smpcats2[j] in list(metaPlot_i[var2]):
                    metaPlot_ij = metaPlot_i[metaPlot_i[var2]==smpcats2[j]]
                    xlist = metaPlot_ij[xn]
                    ylist = metaPlot_ij[yn]
                    ax.scatter(xlist, ylist, label=None, color=colorlist[i], marker=markerlist[jcounter], s=170)

                    if jcounter not in shapeTracker:
                        linesShape[0].append(ax.scatter([], [], label=str(smpcats2[j]), color='black', marker=markerlist[jcounter]))
                        linesShape[1].append(smpcats2[jcounter])
                        shapeTracker.append(jcounter)
                jcounter += 1

            # Here set both legends for color and marker
            ax.legend(linesColor[0], linesColor[1], ncol=1, bbox_to_anchor=(1, 1), title='Sample type', frameon=False, markerscale=1.9, fontsize=18, loc=2)
            from matplotlib.legend import Legend
            leg = Legend(ax, linesShape[0], linesShape[1], bbox_to_anchor=(1, 0.4), title='Electron donor', frameon=False, markerscale=1.9, fontsize=18, loc=2)
            ax.add_artist(leg)

        else: #If there is no var2, change both color and marker with each category in var1
            linesColor[0].append(ax.scatter([], [], label=str(smpcats1[i]), color=colorlist[i], marker=markerlist[i]))
            linesColor[1].append(smpcats1[i])

            xlist = metaPlot_i[xn]
            ylist = metaPlot_i[yn]
            ax.scatter(xlist, ylist, label=None, color=colorlist[i], marker=markerlist[i], s=100)
            ax.legend(linesColor[0], linesColor[1], ncol=1, bbox_to_anchor=(1, 1), title='', frameon=False, markerscale=1.5, fontsize=18)

    ##Put tags at each point
    if tag != 'None':
        for ix in metaPlot.index:
            tagtext = metaPlot.loc[ix, tag]
            tagx = metaPlot.loc[ix, xn]
            tagy = metaPlot.loc[ix, yn]
            ax.annotate(tagtext, (tagx, tagy))

    ax.set_xlabel(xn)
    ax.set_ylabel(yn)
    plt.title(title)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    if savename != 'None':
        plt.savefig(savename)
    plt.show()