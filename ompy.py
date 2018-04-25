import numpy as np
import math
import random
import pickle
import Levenshtein as Lv
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
import ecopy
import pandas as pd
import time
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import seaborn as sns

pd.options.mode.chained_assignment = None  # default='warn'

# Default colors used in plotting
colorlist1 = ['black', 'cyan', 'blue', 'blue', 'blue', 'red', 'darkred', 'darkred', 'darkred', 'grey',
             'darkgrey', 'darkgrey', 'darkgrey', 'gold', 'orange', 'orange', 'orange', 'darkmagenta', 'magenta',
              'fuchsia', 'pink', 'Burlywood', 'green', 'chartreuse', 'Burlywood']

colorlistX = ['blue', 'red', 'gold', 'magenta', 'darkgrey', 'darkorange', 'tomato', 'khaki', 'darkblue', 'brown',
             'hotpink', 'greenyellow', 'green', 'chartreuse', 'Burlywood', 'green', 'chartreuse', 'Burlywood',
              'green', 'chartreuse', 'Burlywood', 'green', 'chartreuse', 'Burlywood']
colorlist2 = ['darkblue', 'dodgerblue', 'crimson', 'red', 'darkgrey', 'lightgrey', 'darkorange', 'gold', 'magenta', 'pink', 'darkgreen', 'lawngreen',
             'darkblue', 'dodgerblue', 'crimson', 'red', 'darkgrey', 'lightgrey', 'darkorange', 'gold', 'magenta',
             'pink', 'darkgreen', 'lawngreen',
             'darkblue', 'dodgerblue', 'crimson', 'red', 'darkgrey', 'lightgrey', 'darkorange', 'gold', 'magenta',
             'pink', 'darkgreen', 'lawngreen']

mycmap = colors.ListedColormap(colorlist2, N=12)
markerlist = ['X', '.', 'v', 's', '*', '.', 'v', 's', '*', '.', 'v', 's', '*', '.', 'v', 's', '*', '.', 'v', 's', '*',
              '.', 'v', 's', 'P', 'D', '*', 'X', 'p', '<', '.', ',', '^', '>', '1', '2', '3', '4', '8', 'h', 'H', '+',
              'x', 'd', '_',
              'o', 'v', 's', 'P', 'D', '*', 'X', 'p', '<', '.', ',', '^', '>', '1', '2', '3', '4', '8', 'h', 'H', '+',
              'x', 'd', '_',]

# ------------------------------------------
# FUNCTIONS FOR LOADING AND SAVING DATA FILES

# Count table, SV names should be in first column, the samples with counts, then taxa starting with Kingdom or Domain
# Return dictionary object with 'tab', 'ra', 'tax', 'seq', and 'meta'. All as pandas dataframes.
def loadFiles(path='None', tab='None', fasta='None', meta='None', sep=','):  # Import file and convert them to suitable format
    if path=='None' or tab=='None':
        print('Error: No path or tab specified')
        return 0

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


    if taxaavailable == 1:
        taxtab = readtab.iloc[:, taxpos:]
        taxtab = taxtab.sort_index()

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
        tabSVs = list(ctab.index)
        corrfastalist = []
        for i in range(len(fastalist)):
            if fastalist[i][0] in tabSVs:
                corrfastalist.append(fastalist[i])
        seqtab = pd.DataFrame(corrfastalist, columns=['SV', 'seq'])
        seqtab = seqtab.set_index('SV')
        seqtab = seqtab.sort_index()

    # Read meta data
    if meta != 'None':
        readmeta = pd.read_csv(path + meta, sep=sep, header=0, index_col=0)

    # Return dictionary object with all dataframes
    out = {}
    out['tab'] = ctab
    out['ra'] = ratab
    if taxaavailable == 1:
        out['tax'] = taxtab
    else:
        print('No taxa found')
    if fasta != 'None':
        out['seq'] = seqtab
    else:
        print('No fasta file found')
    if meta != 'None':
        out['meta'] = readmeta
    else:
        print('No metadata found')
    return out


# Function the takes dictionary object with tab, tax, meta, and seq
# and saves count_tax, meta, and seq fasta files as output to the desired path folder
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


# Function that prints some stats on the data
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
# Subset samples based on metadata
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
    if keep0 == False:
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


# Funtion that subset dataframes based on list of SVs to keep
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


# Function that sets values below a certain relative abundance threshold to zero
# Cutoff specifies the RA cutoff in %
# If ra = 'sample', it is based on the relative abundance in a sample, otherwise it based on SV
def removeLowRA(obj, cutoff=0.1, ra='sample', keep0=False):
    if 'tab' not in obj:
        print('Error removeLowRA, tab not in object')
        return 0
    tab = obj['tab'].copy()
    smp_ra = 100*tab/tab.sum()
    sv_ra = 100*tab.div(tab.sum(axis=1), axis=0)
    if ra == 'sample':
        tab[smp_ra < cutoff] = 0
    else:
        tab[sv_ra < cutoff] = 0

    outobj = {}
    if keep0 != True:
        tab = tab[tab.sum(axis=1)>0]
        outobj['tab'] = tab

        if 'seq' in obj:
            seq = obj['seq']
            seq = seq.loc[tab.index, :]
            outobj['seq'] = seq
        if 'tax' in obj:
            tax = obj['tax']
            tax = tax.loc[tab.index, :]
            outobj['tax'] = tax

    else:
        outobj['tab'] = tab
        if 'seq' in obj:
            outobj['seq'] = obj['seq']
        if 'tax' in obj:
            outobj['tax'] = obj['tax']

    if 'meta' in obj:
        outobj['meta'] = obj['meta']

    outobj['ra'] = 100*tab/tab.sum()

    return outobj


# Function that modifies dataframes based on the number of reads
# All elements with read count <= to cutoff are changed to 0
def removeLowReads(obj, cutoff=1):
    tab = obj['tab']
    ra = obj['ra']
    ra[tab <= cutoff] = 0
    tab[tab <= cutoff] = 0
    tab['sum'] = tab.sum(axis=1)

    out = {}
    tab2 = tab[tab['sum'] > 0]
    ra2 = ra[tab['sum'] > 0]
    out['ra'] = ra2
    if 'tax' in obj.keys():
        tax = obj['tax']
        tax2 = tax[tab['sum'] > 0]
        out['tax'] = tax2
    if 'seq' in obj.keys():
        seq = obj['seq']
        seq2 = seq[tab['sum'] > 0]
        out['seq'] = seq2
    tab2 = tab2.drop('sum', axis=1)
    out['tab'] = tab2
    if 'meta' in obj:
        out['meta'] = obj['meta']
    return out


# Function that subsets SVs based on the fraction of samples in which they are observed
# Cutoff is in % of samples
# If abund=True the most abundant SVs are kept, otherwise the rare ones are kept
def subsetFreq(obj, cutoff=50, abund=True):
    ra = obj['ra']
    tab = obj['tab']
    tax = obj['tax']
    seq = obj['seq']

    tab['freq'] = 100 * (tab.astype(bool).sum(axis=1)) / len(tab.columns)

    if abund == True:  # Remove SV which don't meet freq criteria
        tab2 = tab[tab['freq'] >= cutoff]
        ra2 = ra[tab['freq'] >= cutoff]
        tax2 = tax[tab['freq'] >= cutoff]
        seq2 = seq[tab['freq'] >= cutoff]
        tab2 = tab2.drop('freq', axis=1)
    else:
        tab2 = tab[tab['freq'] < cutoff]
        ra2 = ra[tab['freq'] < cutoff]
        tax2 = tax[tab['freq'] < cutoff]
        seq2 = seq[tab['freq'] < cutoff]
        tab2 = tab2.drop('freq', axis=1)

    out = {}
    out['tab'] = tab2
    out['ra'] = ra2
    out['tax'] = tax2
    out['seq'] = seq2
    if 'meta' in obj:
        out['meta'] = obj['meta']
    return out


# Merges samples based on variable and list from meta data
def mergeSamples(obj, var, slist, keep0=False):
    if 'meta' not in obj.keys():
        print('Metadata missing')
        return 0

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


# Table with counts rarefied based on counts or cutoff percentage
# If method=rabund --> for each sample, check if SV rel abund is > cutoff%, if not - set to 0
# If method=round --> First round, then add or subtract counts if total counts in a sample is not equal to depth
# If method=random --> Randomly pick depth counts of present SVs. Chance of being picked is proportional to number of counts.
def rarefy(tab, method='round', depth='min', cutoff=0.1):
    tab = tab.applymap(int) #Make sure table elements are integers
    if depth == 'min':  # Define read depth
        reads = min(tab.sum())
    else:
        reads = depth

    if method == 'round':
        samples = tab.columns.tolist()
        svs = tab.index.tolist()
        rtab = tab.copy()
        for smp in samples:
            smpsum = sum(tab[smp])
            if smpsum < reads:  # Remove sample if sum of reads less than read depth
                rtab = rtab.drop(smp, axis=1)
                continue
            frac = tab.loc[:, smp] * reads / smpsum
            avrundad = frac.apply(round)
            diffs = frac - avrundad  # Round calculatad read numbers
            if sum(avrundad) < reads:  # Correct if too few
                addNR = int(reads - sum(avrundad))
                diffs = diffs[diffs > 0]
                addSV = random.sample(list(diffs.index), addNR)
                avrundad[addSV] = avrundad[addSV] + 1
                rtab[smp] = avrundad
            elif sum(avrundad) > reads:  # Correct if too many
                addNR = int(sum(avrundad) - reads)
                diffs = diffs[diffs < 0]
                addSV = random.sample(list(diffs.index), addNR)
                avrundad[addSV] = avrundad[addSV] - 1
                rtab[smp] = avrundad
            else:
                rtab[smp] = avrundad
        return rtab

    elif method == 'random':
        rftab = pd.DataFrame(0, index=tab.index, columns=tab.columns)
        for smp in tab.columns:
            if tab[smp].sum() < reads:
                rftab = rftab.drop([smp], axis=1)
                #print(smp, ' dropped during rarefy')
                continue
            else:
                temparr = tab[smp][tab[smp] > 0]
                longlist = []
                for sv in temparr.index:
                    longlist = longlist + [sv] * int(temparr[sv])
                randsample = random.sample(longlist, reads)
                uniquelist = np.unique(randsample, return_counts=True)
                rftab.loc[uniquelist[0], smp] = uniquelist[1]
        return rftab

    elif method == 'rabund':
        rtab = tab.copy()
        for smp in tab.columns:
            ra = 100 * rtab[smp] / rtab[smp].sum()
            rtab[smp][ra < cutoff] = 0
        return rtab

    else:
        print('Error in Rarefy input')

# Table with rarefied counts, this method is much faster than rarefy with random method above
# Randomly pick depth counts of present SVs. Chance of being picked is proportional to number of counts (i.e. with replacement).
def rarefy2(tab, depth='min', seed=0):
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
        if tab[s].sum()<depth:
            rftab = rftab.drop(s, axis=1)
            continue
        else:
            p = tab[s]/tab[s].sum()
            choice = prng.choice(nvar, depth, p=p)
            rftab[s] = np.bincount(choice, minlength=nvar)
    return rftab

# Table with rarefied counts, this method is much faster than rarefy with random method above
# Randomly pick depth counts of present SVs. This one is without replacement.
def rarefy3(tab, depth='min', seed=42):
    steps = 20
    prng = np.random.RandomState(seed)  # reproducible results
    noccur = tab.sum()
    nvar = len(tab.index)  # number of SVs

    ## Set read depth
    if depth == 'min':
        depth = int(np.min(noccur))
    else:
        depth = int(depth)

    rftab = tab.copy()
    for s in tab.columns:  # for each sample
        valuelist = tab[s]
        if valuelist.sum() < depth:
            rftab = rftab.drop(s, axis=1)
            continue
        else:
            rftab[s] = 0
            for i in range(steps):
                p = valuelist / valuelist.sum()
                choice = prng.choice(nvar, int(depth/steps), p=p)
                binned_choice = np.bincount(choice, minlength=nvar)
                rftab[s] = rftab[s] + binned_choice
                valuelist = valuelist - binned_choice
                valuelist[valuelist<0] = 0
    return rftab

# Group samples in a results series based on meta variable
# Take averages from each group
def sampleSeriesAverages(ser, meta, var, slist):
    slist = sorted(list(set(slist)))
    avglist = []
    for s in slist:
        smps = meta[meta[var]==s].index.tolist()
        avglist.append(ser.loc[smps].mean())
    return pd.Series(avglist, index=slist)

# ------------------------------------------
# FUNCTIONS FOR VISUALISING TAXA IN HEATMAP

# Group SVs based on taxa
# Returns only tab and ra
def groupbyTaxa(obj, level='Genus', extra='none'):
    # Clean up tax data frame
    tax = obj['tax']
    tax = tax.fillna(0)

    # Get index of levels and extra
    taxanameslist = tax.columns.tolist() #List with Kingdom, Phylum etc.
    levelix = taxanameslist.index(level)
    extraix = taxanameslist.index(extra)

    # Change Otu into SV
    indexlist = list(tax.index)
    if indexlist[0][:3] == 'Otu': #Check if Otu are used in index
        newnames = []
        for i in range(len(indexlist)):
            newnames.append('SV' + indexlist[i][3:])
        indexlist = newnames

    # Put the SV name in all empty spots in dataframe called taxSV
    taxSV = tax.copy()
    for s in range(len(taxSV.columns)):
        for i in range(len(taxSV.index)):
            if taxSV.iloc[i, s] == 0:
                taxSV.iloc[i, s] = indexlist[i]
    # Change all 0 in tax to lowest determined taxa level
    taxNames = tax.copy()
    for s_nr in range(1, len(taxanameslist)):
        s1 = taxanameslist[s_nr]; s0 = taxanameslist[s_nr-1]
        taxNames[s1][tax[s1] == 0] = taxNames[s0][tax[s1] == 0]

    # Create names to use in output
    if extra in tax.columns and level in tax.columns:
        taxSV['Name'] = taxSV[extra] + '; ' + taxSV[level]
        for i in range(len(tax.index)): #Go through names and check so none is SV-SV
            if tax[extra].iloc[i] == 0:
                taxSV['Name'].iloc[i] = str(taxNames[extra].iloc[i])+'; '+ str(taxSV[level].iloc[i])
    elif level in tax.columns:
        taxSV['Name'] = taxSV[level]
        for i in range(len(tax.index)): #Go through names and check so none is only SV
            if tax[level].iloc[i] == 0:
                taxSV['Name'].iloc[i] = taxNames[level].iloc[i]+'; '+taxSV[level].iloc[i]
    else:
        print('Error in GroupbyTaxa, level is not in tax table')
        return 0

    tab = obj['tab']
    tab['Name'] = list(taxSV['Name'])
    ra = obj['ra']
    ra['Name'] = list(taxSV['Name'])

    tab = tab.set_index(['Name'])
    ra = ra.set_index('Name')
    tab = tab.groupby(tab.index).sum()
    ra = ra.groupby(ra.index).sum()

    out = {}
    out['tab'] = tab
    out['ra'] = ra
    return out


# Returns a dataframe with rounded relative abundance values for heatmap
def heatmapRAdataLabels(ra):
    for r in ra.index:
        for c in ra.columns:
            value = float(ra.loc[r, c])
            if value < 0.1 and value > 0:
                ra.loc[r, c] = -0.1
            elif value < 10 and value >= 0.1:
                ra.loc[r, c] = round(value, 1)
            elif value > 99:
                ra.loc[r, c] = int(99)
            elif value >= 10:
                ra.loc[r, c] = round(value)
            else:
                ra.loc[r, c] = 0
    return ra


# For each sample retain top number SVs
def listTopSVs(obj, number=20):
    tab = obj['tab'].copy()
    ra = obj['ra'].copy()
    cols = tab.columns.tolist()
    retain = []
    # This goes through the samples and picks out the top SVs in each sample
    for c in cols:
        tab = tab.sort_values(by=[c], ascending=False)
        svs = tab.index.tolist()[:number]
        retain = retain + svs
    # This identifies top number SVs based on max RA in a sample
    ra['max'] = ra.max(axis=1)
    ra = ra.sort_values(by=['max'], ascending=False)
    svs = ra.index.tolist()[:number]
    retain = retain + svs
    # Get the uniqe SVs in retain and return as list
    retain = set(retain)
    return list(retain)

# Plots heatmap, order is the heading in metadata for the column that specifies logical order of samples
def plotHeatmap(obj, var, taxalevels=['Phylum', 'Genus'], order='None', SVs2plot=20, savename='None'):
    #Merge samples based on var
    merged_obj = mergeSamples(obj, var=var, slist=obj['meta'].loc[:, var])

    ## Make sure samples are in the right order in meta data
    if order != 'None':
        md = merged_obj['meta']
        md[order] = md[order].astype(float)
        md = md.sort_values(by=order)
        logiclist = []
        [logiclist.append(item) for item in md[var] if item not in logiclist]
        merged_obj['meta'] = md

    ## Groupby taxa
    taxa_obj = groupbyTaxa(merged_obj, level=taxalevels[1], extra=taxalevels[0])

    # Subset for top taxa
    topSVs = listTopSVs(taxa_obj, number=SVs2plot)
    retained = subsetSVs(taxa_obj, topSVs)

    # Print heatmap
    table = retained['ra']
    if order != 'None':
        table = table.loc[:, logiclist]

    table['max'] = table.max(axis=1)
    table = table.sort_values(by=['max'], ascending=False)
    table = table.iloc[:SVs2plot, :]
    table = table.drop(['max'], axis=1)

    table['avg'] = table.mean(axis=1)
    table = table.sort_values(by=['avg'], ascending=True)
    table = table.drop(['avg'], axis=1)

    #Fix datalabels
    datalab = heatmapRAdataLabels(table)

    #Plot
    fig, ax = plt.subplots(figsize=(14, 10))
    sns.set(font_scale=1.2)
    sns.heatmap(table, annot=datalab, cmap='Reds', linewidths=0.5, robust=True, cbar=False, ax=ax)
    plt.xticks(rotation=90)
    plt.setp(ax.get_xticklabels(), fontsize=16)
    ax.set_ylabel('')
    plt.setp(ax.get_yticklabels(), fontsize=16)
    plt.tight_layout()
    if savename != 'None':
        plt.savefig(savename)
    plt.show()


    # plt.figure(figsize=(14, 10))
    # sns.set(font_scale=1.1)
    # sns.heatmap(table, annot=datalab, cmap='Reds', linewidths=0.5, robust=True, cbar=False, ax=ax)
    # plt.xticks(rotation=90)
    # plt.ylabel('')
    # plt.setp(ax.get_yticklabels(), fontsize=14)
    # plt.tight_layout()
    # if savename != 'None':
    #     plt.savefig(savename)
    # plt.show()

#blablabl
# ------------------------------------------
# ALPHA AND BETA DIVERSITY FUNCTIONS

# Calculates Hill numbers of order q based on Jost 2006
def naiveDivAlpha(tab, q=1):
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


# Function to calculate pairwise beta diversity matrix
# Based on Jost 2007, Ecology 88(10), p2427
def naiveDivBeta(tab, q=1):
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
    return outdf

# Calculate Bray-Curtis dissimilarity matrix
def bray(tab):
    ra = tab / tab.sum()
    smplist = list(tab.columns)
    outdf = pd.DataFrame(0, index=smplist, columns=smplist)
    for smp1nr in range(len(smplist) - 1):
        smp1 = smplist[smp1nr]
        for smp2nr in range(smp1nr + 1, len(smplist)):
            smp2 = smplist[smp2nr]
            brayvalue = 1-(ra.loc[:, [smp1, smp2]].min(axis=1).sum())
            outdf.loc[smp1, smp2] = brayvalue; outdf.loc[smp2, smp1] = brayvalue
    return outdf

# Calculates Pielou or Simpson evenness index
def evenness(tab, index='Pielou'):  # Dataframe with sample names as row index and evenness parameters as values
    if index == 'Pielou':
        Rtab = tab.copy()
        R = naiveDivAlpha(Rtab, q=0)
        Hmax = R.apply(math.log)

        ra = tab / tab.sum()
        mask = ra != 0
        raLn = ra[mask] * np.log(ra[mask])
        H = -1 * raLn.sum()
        return H / Hmax
    elif index == "Simpson":
        ra = tab / tab.sum()
        ra2 = ra.mul(ra)
        return 1 - ra2.sum()


# ------------------------------------------
# PHYLOGENETIC DIVERSITY AND RELATEDNESS FUNCTIONS

# Calculates a matrix for pairwise distances between all SV sequences
def phylDistMat(seq, savename='PhylDistMat.pickle'):
    svnames = list(seq.index)
    df = pd.DataFrame(0, index=svnames, columns=svnames)
    for i in range(len(svnames) - 1):
        if i % 10 == 0:
            print(i)
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
    with open(savename, 'wb') as handle:
        pickle.dump(df, handle)
    print('Finished printing Phyldistmat as pickle file')


# Returns Rao's quadratic entropy, sum(sum(dij*pi*pj))
# Used in Chiu's phylogenetic diversity functions
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


# Calculates phylogenetic (functional) diversity according to Chiu's framework in PLOS One, 2014
def phylDivAlpha(tab, distmat, q=0, index='FD'):
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
    if index == 'D':
        return outdf
    elif index == 'MD':
        return outdf.mul(Qframe)
    elif index == 'FD':
        MD = outdf.mul(Qframe)
        return outdf.mul(MD)


# Calculates phylogenetic (functional) beta diversity according to Chiu's framework in PLOS One, 2014
def phylDivBeta(tab, distmat, q=1, index='FD'):
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
                outD.loc[smp1, smp2] = Dg / Da; outD.loc[smp2, smp1] = Dg / Da

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
                outD.loc[smp1, smp2] = Dg / Da; outD.loc[smp2, smp1] = Dg / Da
    if index == 'FD':
        outFD = outD.pow(2)
        return outFD
    elif index == 'D':
        return outD


# Converts beta value to distances, specify q and type associated with the beta (assuming pairwise)
def beta2Dist(beta, q=1, type='naive'):
    beta = beta.applymap(float)
    dist = beta.copy()
    mask = beta > 0
    if type == 'naive':
        if q == 1:
            dist[mask] = 1-(math.log(2)-np.log(beta[mask]))/math.log(2)
        else:
            dist[mask] = 1-((1/beta[mask])**(q-1)-0.5**(q-1))/(1-0.5**(q-1))

    elif type == 'phyl':
        if q == 1:
            dist[mask] = np.log(beta[mask]) / (2 * math.log(2))
        else:
            dist[mask] = 1-((2**(2*(1-q))-beta[mask].pow(1-q))/(2**(2*(1-q))-1))
    return dist


# Calculates mean nearest neighbour distance for each sample
# Weighs SVs based on abundance, q value determines how much
def aMNND(tab, distmat, q):
    output = pd.Series(0, index=tab.columns)
    for smp in tab.columns:
        nonzeroSVs = tab[smp][tab[smp]>0].index
        tabpow = tab.loc[nonzeroSVs,smp].pow(q)
        ra = tabpow/tabpow.sum()

        nonzeroDM = distmat[distmat > 0]

        output.loc[smp] = nonzeroDM.loc[nonzeroSVs,nonzeroSVs].min(axis=1).mul(ra).sum()
    return output


# Calculated beta mean nearest taxonomic distance as matrix - NEED TO FIX!
def bMNTD(tab, distmat):
    ra = tab / tab.sum()
    smplist = list(tab.columns)

    outdf = pd.DataFrame(0, index=smplist, columns=smplist)
    for s1nr in range(len(smplist) - 1):
        for s2nr in range(s1nr + 1, len(smplist)):
            s1 = smplist[s1nr]
            s2 = smplist[s2nr]
            s1SVs = list(tab[tab[s1] > 0].index)
            s2SVs = list(tab[tab[s2] > 0].index)

            dm = distmat.loc[s1SVs, s2SVs]
            ntd = 0
            ntd = ntd + dm.min(axis=1).mul(ra.loc[s1SVs, s1]).sum()
            ntd = ntd + dm.min(axis=0).mul(ra.loc[s2SVs, s2]).sum()
            ntd = 0.5 * ntd
            outdf.loc[s1, s2] = ntd
            outdf.loc[s2, s1] = ntd
    return outdf


# Plots PCoA for distance matrix
def plotPCoA(dist, meta, var1='None', var2='None', order='None', title='', savename='None'):
    # Do the PCoA in ecopy and get results for PC1 and PC2 into dataframe
    pcs = ecopy.pcoa(dist)
    pc_values = ecopy.pcoa.summary(pcs)
    pc1_perc = round(100 * pc_values.iloc[1, 0], 1);
    xn = 'PC1 (' + str(pc1_perc) + '%)'
    pc2_perc = round(100 * pc_values.iloc[1, 1], 1);
    yn = 'PC2 (' + str(pc2_perc) + '%)'

    koord = pcs.biplot(coords=True)['Objects']
    smplist = dist.index
    pcoadf = pd.DataFrame(koord, index=smplist, columns=[xn, yn])

    # Combine pcoa results with meta data
    meta[xn] = pcoadf[xn]
    meta[yn] = pcoadf[yn]
    meta = meta[meta[xn].notnull()]
    if order != 'None':
        meta = meta.sort_values(by=[order])

    if var1 == 'None' and var2 == 'None':
        return 'Error, no variables in input'
    if var1 != 'None':
        smpcats1 = []
        [smpcats1.append(item) for item in meta[var1] if item not in smpcats1]
    if var2 != 'None':
        smpcats2 = []
        [smpcats2.append(item) for item in meta[var2] if item not in smpcats2]

    # Do the plotting
    metaPlot = meta.loc[:, [xn,yn, var1, var2]]

    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(10, 6))

    linesColor = [[], []]
    linesShape = [[], []]
    shapeTracker = []

    for i in range(len(smpcats1)):
        metaPlot_i = metaPlot[metaPlot[var1] == smpcats1[i]]

        if var2 != 'None':
            linesColor[0].append(ax.scatter([], [], label=str(smpcats1[i]), color=colorlist1[i]))
            linesColor[1].append(smpcats1[i])

            jcounter = 0
            for j in range(len(smpcats2)):
                if smpcats2[j] in list(metaPlot_i[var2]):
                    metaPlot_ij = metaPlot_i[metaPlot_i[var2]==smpcats2[j]]
                    xlist = metaPlot_ij[xn]
                    ylist = metaPlot_ij[yn]
                    ax.scatter(xlist, ylist, label=None, color=colorlist1[i], marker=markerlist[jcounter], s=70)

                    if jcounter not in shapeTracker:
                        linesShape[0].append(ax.scatter([], [], label=str(smpcats2[j]), color='black', marker=markerlist[jcounter]))
                        linesShape[1].append(smpcats2[jcounter])
                        shapeTracker.append(jcounter)
                jcounter += 1

            ax.legend(linesColor[0], linesColor[1], ncol=1, bbox_to_anchor=(0.95, 1.0), title='Enrichment', frameon=False)
            from matplotlib.legend import Legend
            leg = Legend(ax, linesShape[0], linesShape[1], bbox_to_anchor=(1.2, 1.0), title='Day', frameon=False)
            ax.add_artist(leg)

        else:
            linesColor[0].append(ax.scatter([], [], label=str(smpcats1[i]), color=colorlist1[i], marker=markerlist[i]))
            linesColor[1].append(smpcats1[i])

            xlist = metaPlot_i[xn]
            ylist = metaPlot_i[yn]
            ax.scatter(xlist, ylist, label=None, color=colorlist1[i], marker=markerlist[i], s=100)
    lgnd = ax.legend(linesColor[0], linesColor[1], ncol=1, bbox_to_anchor=(1.0, 1.05), title='', frameon=False, markerscale=1.5, fontsize=14)
    # for lgi in range(len(smpcats1)):
    #     lgnd.legendHandles[lgi]._sizes = [120]

    ax.set_xlabel(xn)
    ax.set_ylabel(yn)
    plt.title(title)
    plt.tight_layout(rect=[0, 0, 0.75, 1])
    if savename != 'None':
        plt.savefig(savename)
    plt.show()


# ----------------------------------
# COMPARE TECHNICAL REPLICATES -#

# Plots the mean relative abundance vs the coefficient of variation for SVs
# Create list of pairwise comparisons of RA for all SVs and samples
def plotRepComp(tab, xlims=[0, 100], ylims=[0, 100], savename='plotrepcomp'):
    ctab = tab.copy()
    smps = ctab.columns
    ra = 100 * ctab / ctab.sum()
    ra['mean'] = ra.mean(axis=1)
    ra['std'] = ra[smps].std(axis=1)
    ra['cov'] = 100 * ra['std'] / ra['mean']
    ra['min'] = ra[smps].min(axis=1)
    ra['max'] = ra[smps].max(axis=1)
    ra['not_seen'] = 100*ra[smps][ra==0].count(axis=1)/len(tab.columns)
    ctab['min'] = ctab.min(axis=1)
    ctab['max'] = ctab.max(axis=1)

    plt.rcParams.update({'font.size': 14})
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)  # Create figure
    ax1.plot(ra['max'], ra['not_seen'], '.')
    plt.xlabel('Mean rel. abund. (%)')
    plt.xlabel('Max rel. abund. (%)')
    #plt.ylabel('Coeff. of variation (%)')
    plt.ylabel('Fraction of samples not detected (%)')
    plt.xlim(xlims[0], xlims[1])
    plt.ylim(ylims[0], ylims[1])
    plt.tight_layout()
    #plt.savefig(savename)
    plt.show()

# Plots a rarefaction curves based on random subsampling
# the colors of the lines are based on var and slist from metadata
def plotRarefaction(obj, var, slist='all', step=2000, savename='None'):
    tab = obj['tab']; tab = tab.applymap(int)
    meta = obj['meta']

    if slist == 'all':
        slist = meta[var]

    submeta = meta.loc[meta[var].isin(slist)]
    smplist = list(submeta.index)
    subtab = tab.loc[:, smplist]

    maxcount = int(max(subtab.sum()))
    plottab = [['count']+smplist]
    numberofsteps = int(maxcount/step); counter = 0
    for i in range(step, maxcount, step):
        counter+=1
        if counter%5 == 0:
            print(counter, ' out of ', numberofsteps, ' steps')
        raretab = rarefy3(subtab, depth=i)
        svcount = naiveDivAlpha(raretab, 0)
        templist=[i]
        for s in smplist:
            if s in svcount.index:
                templist = templist+[svcount[s]]
            else:
                templist = templist+['NaN']
        plottab.append(templist)
    plottab = pd.DataFrame(plottab[1:], columns=plottab[0])
    plottab.to_csv(savename+'.csv')

    #Plot
    colorcheck = []
    plt.rcParams.update({'font.size': 20})
    plt.figure(figsize=(10, 6))
    for s in smplist:
        stype = submeta.loc[s, var]
        if stype not in colorcheck:
            colorcheck.append(stype)
            plt.plot(plottab['count'],plottab[s], label=stype, color=colorlist[len(colorcheck)-1])
        else:
            plt.plot(plottab['count'], plottab[s], label='_nolegend_', color=colorlist[len(colorcheck) - 1])
    plt.legend(loc='best')
    plt.xlabel('Number of reads')
    plt.ylabel('Richness')
    plt.tight_layout()
    if savename != 'None':
        plt.savefig(savename+'.png')
    plt.show()

def plotRankAbundance(tab, slist):
    fig, ax = plt.subplots(1, figsize=(8,5))
    plt.rcParams.update({'font.size':14})

    for s in slist:
        ra = 100*tab[s]/tab[s].sum()
        ra = ra[ra>0]
        ra = ra.sort_values(ascending=False)
        xaxis = np.arange(1,len(ra)+1)
        ax.plot(xaxis, ra, marker='.', lw=0, label=s)
    ax.set_xlim(0,10)
    ax.legend()
    plt.tight_layout()
    plt.show()

# ------------------------------------
# FUNCTIONS FOR CORE COMMUNITY

# Returns list of SV present with rel abund above cutoff % in all samples in slist
def commonSVs(obj, var, slist, cutoff=0.05):  # List of common SVs
    ra = obj['ra'];
    md = obj['meta']
    SVdict = {}
    for smptype in slist:
        samples = md.index[md[var] == smptype].tolist()
        ra['min'] = ra[samples].min(axis=1)
        SVsinsmp = ra.index[ra['min'] > cutoff].tolist()
        SVdict[smptype] = SVsinsmp
    if len(slist) > 1:
        for i in range(0, len(slist) - 1):
            out = list(set(SVdict[slist[i]]).intersection(SVdict[slist[i + 1]]))
            SVdict[slist[i + 1]] = out
        return out
    else:
        print('More items needed in slist -- Common SVs')


# Finds SV present in only that sample type (over cutoff2) and absent in other sample types (below cutoff2)
# Returns dictionary with list for each sample type
def uniqueSVs(obj, var, slist, cutoff1=0.2, cutoff2=0.01):
    ra = obj['ra'];
    md = obj['meta']
    SVpres = {};
    SVabs = {}
    for smptype in slist:
        samples = md.index[md[var] == smptype].tolist()
        ra['min'] = ra[samples].mean(axis=1)
        SVsinsmp = ra.index[ra['min'] > cutoff1].tolist()
        ra['max'] = ra[samples].mean(axis=1)
        SVsnotin = ra.index[ra['max'] <= cutoff2].tolist()
        SVpres[smptype] = SVsinsmp;
        SVabs[smptype] = SVsnotin
    if len(slist) > 1:
        for i in range(len(slist)):
            for j in range(len(slist)):
                if i == j:
                    continue
                else:
                    out = list(set(SVpres[slist[i]]).intersection(SVabs[slist[j]]))
                    SVpres[slist[i]] = out
        return SVpres
    else:
        print('More items needed in slist -- Common SVs')


# -------------------------
# NULL MODELS

def neutralAnimation(path):

    totalpool = 625; side1 = 25; side2 = 25
    exchange = 220

    species = pd.DataFrame(1, index=range(1,11), columns=['R1','R2'])
    randomspecies = species.sample(totalpool, replace=True)
    R1pool = pd.Series(randomspecies.index)
    R2pool = pd.Series(randomspecies.index)

    xtime = []
    R1D0 = []; R2D0 = []; R1D1 = []; R2D1 = []
    R12beta0 = []; R12beta1 = []

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10,10))


    def make_frame(i):
        nonlocal R1pool, R2pool, species
        nonlocal xtime, R1D0, R2D0, R1D1, R2D1, R12beta0, R12beta1

        #Identify index of those to die
        toDie1 = list(R1pool.sample(exchange, replace=False).index)
        toDie2 = list(R2pool.sample(exchange, replace=False).index)

        #Set those index to 0 and return unique counts of remaining
        R1pool[toDie1] = 0; R2pool[toDie2] = 0
        R1uni = np.unique(R1pool, return_counts=True)
        R2uni = np.unique(R2pool, return_counts=True)

        #Add those unique counts to species abundance
        species.loc[:,['R1','R2']] = 0
        species.loc[R1uni[0][1:], 'R1'] = R1uni[1][1:]
        species.loc[R2uni[0][1:], 'R2'] = R2uni[1][1:]

        #Get random species (of remaining) to replace those that died
        randomspecies1 = species.sample(exchange, replace=True, weights=species['R1'])
        randomspecies2 = species.sample(exchange, replace=True, weights=species['R2'])

        #Add those random to pool
        R1pool[toDie1] = randomspecies1.index
        R2pool[toDie2] = randomspecies2.index

        #Make matrix
        mat1 = R1pool.values.reshape(side1, side2)
        mat2 = R2pool.values.reshape(side1, side2)

        #Update species df again
        R1uni = np.unique(R1pool, return_counts=True)
        R2uni = np.unique(R2pool, return_counts=True)
        species.loc[:,['R1','R2']] = 0
        species.loc[R1uni[0], 'R1'] = R1uni[1]
        species.loc[R2uni[0], 'R2'] = R2uni[1]

        #Calculate diversity values
        xtime.append(i)
        alpha0 = naiveDivAlpha(species, 0)
        alpha1 = naiveDivAlpha(species, 1)
        R1D0.append(alpha0['R1']); R2D0.append(alpha0['R2'])
        R1D1.append(alpha1['R1']); R2D1.append(alpha1['R2'])

        beta0 = naiveDivBeta(species, 0)
        beta1 = naiveDivBeta(species, 1)
        dis0 = beta2Dist(beta0, 0, 'naive')
        dis1 = beta2Dist(beta1, 1, 'naive')
        R12beta0.append(dis0.loc['R1', 'R2'])
        R12beta1.append(dis1.loc['R1', 'R2'])

        ax[0,0].clear()
        ax[0,1].clear()
        ax[1,0].clear()
        ax[1,1].clear()

        ax[0,0].pcolor(mat1, cmap=mycmap, vmin=0, vmax=11)
        ax[0,1].pcolor(mat2, cmap=mycmap, vmin=0, vmax=11)
        ax[0,0].set_title('R1')
        ax[0,1].set_title('R2')

        ax[1,0].plot(xtime, R1D0, color='blue', lw=1, label='R1_D0')
        ax[1,0].plot(xtime, R2D0, color='green', lw=0.8, label='R2_D0')
        ax[1,0].plot(xtime, R1D1, color='red', lw=1, label='R1_D1')
        ax[1,0].plot(xtime, R2D1, color='gold', lw=0.8, label='R2_D1')
        ax[1,1].plot(xtime, R12beta0, color='blue', lw=1, label='dis0')
        ax[1,1].plot(xtime, R12beta1, color='red', lw=0.8, label='dis1')
        ax[1,0].set_ylabel('Diversity')
        ax[1,0].set_ylim(0,12)
        ax[1,1].set_ylabel('Dissimilarity')
        ax[1,1].set_ylim(0,1)
        ax[1,0].legend(loc=3)
        ax[1,1].legend(loc=2)
        plt.tight_layout()
        #plt.pause(0.1)

        return mplfig_to_npimage(fig)

    animation = VideoClip(make_frame, duration=30)
    animation.write_gif(path+'Neutral_animation.gif', fps=5)

def neutralModel(tab, startsample, iters=1000):
    samplecount = tab[startsample].sum()

    output = [['iter', 'dt_0ND', 'dt_1ND', 'dt_2ND', 'dx_0ND', 'dx_1ND', 'dx_2ND', 'a_0ND', 'a_1ND', 'a_2ND']]
    df = pd.DataFrame({'c1': tab[startsample], 'c2': tab[startsample]})

    prevc1 = df['c1']
    prevc2 = df['c2']
    dtemp = df.copy()
    dtemp['prevc1'] = prevc1
    dtemp['prevc2'] = prevc2

    alfa0 = naiveDivAlpha(dtemp, 0)
    alfa1 = naiveDivAlpha(dtemp, 1)
    alfa2 = naiveDivAlpha(dtemp, 2)
    beta0 = naiveDivBeta(dtemp, 0)
    beta1 = naiveDivBeta(dtemp, 1)
    beta2 = naiveDivBeta(dtemp, 2)
    dist0 = beta2Dist(beta0, 0, 'naive')
    dist1 = beta2Dist(beta1, 1, 'naive')
    dist2 = beta2Dist(beta2, 2, 'naive')

    outi = 0
    dt0 = np.average([dist0.loc['c1', 'prevc1'], dist0.loc['c2', 'prevc2']])
    dt1 = np.average([dist1.loc['c1', 'prevc1'], dist1.loc['c2', 'prevc2']])
    dt2 = np.average([dist2.loc['c1', 'prevc1'], dist2.loc['c2', 'prevc2']])
    dx0 = dist0.loc['c1', 'c2']
    dx1 = dist1.loc['c1', 'c2']
    dx2 = dist2.loc['c1', 'c2']
    a0 = np.average(alfa0.values)
    a1 = np.average(alfa1.values)
    a2 = np.average(alfa2.values)
    output.append([outi, dt0, dt1, dt2, dx0, dx1, dx2, a0, a1, a2])

    df[df == 0] = 0.5
    ra = df / df.sum()

    for i in range(1, iters + 1):
        randnrs1 = np.random.rand(len(ra.index)) + 0.5
        randnrs2 = np.random.rand(len(ra.index)) + 0.5
        ra['c1'] = ra['c1'] * randnrs1
        ra['c2'] = ra['c2'] * randnrs2
        ra['c1'] = ra['c1'] / ra['c1'].sum()
        ra['c2'] = ra['c2'] / ra['c2'].sum()

        df['c1'] = ra['c1'] * samplecount
        df['c2'] = ra['c2'] * samplecount
        df = df.applymap(round)
        for c in ['c1', 'c2']:
            diff = int(samplecount - df[c].sum())
            if diff > 0:
                indx = list(df[c].sample(n=diff, weights=ra[c]).index)
                df.loc[indx, c] = df.loc[indx, c] + 1
            elif diff < 0:
                indx = list(df[df[c] > 0].index)
                np.random.shuffle(indx)
                df.loc[indx[:-diff], c] = df.loc[indx[:-diff], c] - 1

        dtemp = df.copy()
        dtemp['prevc1'] = prevc1
        dtemp['prevc2'] = prevc2

        alfa0 = naiveDivAlpha(dtemp, 0)
        alfa1 = naiveDivAlpha(dtemp, 1)
        alfa2 = naiveDivAlpha(dtemp, 2)
        beta0 = naiveDivBeta(dtemp, 0)
        beta1 = naiveDivBeta(dtemp, 1)
        beta2 = naiveDivBeta(dtemp, 2)
        dist0 = beta2Dist(beta0, 0, 'naive')
        dist1 = beta2Dist(beta1, 1, 'naive')
        dist2 = beta2Dist(beta1, 2, 'naive')

        outi = i
        dt0 = np.average([dist0.loc['c1', 'prevc1'], dist0.loc['c2', 'prevc2']])
        dt1 = np.average([dist1.loc['c1', 'prevc1'], dist1.loc['c2', 'prevc2']])
        dt2 = np.average([dist2.loc['c1', 'prevc1'], dist2.loc['c2', 'prevc2']])
        dx0 = dist0.loc['c1', 'c2']
        dx1 = dist1.loc['c1', 'c2']
        dx2 = dist2.loc['c1', 'c2']
        a0 = np.average(alfa0.values)
        a1 = np.average(alfa1.values)
        a2 = np.average(alfa2.values)
        output.append([outi, dt0, dt1, dt2, dx0, dx1, dx2, a0, a1, a2])
        prevc1 = df['c1']
        prevc2 = df['c2']

        if i % 50 == 0:
            print(i, '  ', diff)

    outdf = pd.DataFrame(output[1:], columns=output[0])
    outdf.to_csv('Neutralmodel.csv')

    plt.figure(figsize=(10, 8))
    plt.plot(outdf['iter'], outdf['dt_0ND'], label='dt_0ND')
    plt.plot(outdf['iter'], outdf['dt_1ND'], label='dt_1ND')
    plt.plot(outdf['iter'], outdf['dt_2ND'], label='dt_2ND')
    plt.legend(loc='best')
    plt.show()

    plt.figure(figsize=(10, 8))
    plt.plot(outdf['iter'], outdf['dx_0ND'], label='dx_0ND')
    plt.plot(outdf['iter'], outdf['dx_1ND'], label='dx_1ND')
    plt.plot(outdf['iter'], outdf['dx_2ND'], label='dx_2ND')
    plt.legend(loc='best')
    plt.show()

    plt.figure(figsize=(10, 8))
    plt.plot(outdf['iter'], outdf['a_0ND'], label='a_0ND')
    plt.plot(outdf['iter'], outdf['a_1ND'], label='a_1ND')
    plt.plot(outdf['iter'], outdf['a_2ND'], label='a_2ND')
    plt.legend(loc='best')
    plt.show()


def deterministicModel(tab, startsample, iters=3000):
    samplecount = tab[startsample].sum()

    output = [['iter', 'dt_0ND', 'dt_1ND', 'dt_2ND', 'dx_0ND', 'dx_1ND', 'dx_2ND', 'a_0ND', 'a_1ND', 'a_2ND']]
    df = pd.DataFrame({'c1': tab[startsample], 'c2': tab[startsample]})

    prevc1 = df['c1']
    prevc2 = df['c2']
    dtemp = df.copy()
    dtemp['prevc1'] = prevc1
    dtemp['prevc2'] = prevc2

    alfa0 = naiveDivAlpha(dtemp, 0)
    alfa1 = naiveDivAlpha(dtemp, 1)
    alfa2 = naiveDivAlpha(dtemp, 2)
    beta0 = naiveDivBeta(dtemp, 0)
    beta1 = naiveDivBeta(dtemp, 1)
    beta2 = naiveDivBeta(dtemp, 2)
    dist0 = beta2Dist(beta0, 0, 'naive')
    dist1 = beta2Dist(beta1, 1, 'naive')
    dist2 = beta2Dist(beta2, 2, 'naive')

    outi = 0
    dt0 = np.average([dist0.loc['c1', 'prevc1'], dist0.loc['c2', 'prevc2']])
    dt1 = np.average([dist1.loc['c1', 'prevc1'], dist1.loc['c2', 'prevc2']])
    dt2 = np.average([dist2.loc['c1', 'prevc1'], dist2.loc['c2', 'prevc2']])
    dx0 = dist0.loc['c1', 'c2']
    dx1 = dist1.loc['c1', 'c2']
    dx2 = dist2.loc['c1', 'c2']
    a0 = np.average(alfa0.values)
    a1 = np.average(alfa1.values)
    a2 = np.average(alfa2.values)
    output.append([outi, dt0, dt1, dt2, dx0, dx1, dx2, a0, a1, a2])

    df[df == 0] = 0.5
    ra = df / df.sum()

    randnrs1 = np.random.rand(len(ra.index)) + 0.5
    randnrs2 = np.random.rand(len(ra.index)) + 0.5
    randnrs3 = np.random.rand(len(ra.index)) + 0.5

    for i in range(1, iters + 1):
        if i < iters / 3:
            ra['c1'] = ra['c1'] * randnrs1
            ra['c2'] = ra['c2'] * randnrs1
            ra['c1'] = ra['c1'] / ra['c1'].sum()
            ra['c2'] = ra['c2'] / ra['c2'].sum()
        elif i < 2 * iters / 3:
            ra['c1'] = ra['c1'] * randnrs2
            ra['c2'] = ra['c2'] * randnrs2
            ra['c1'] = ra['c1'] / ra['c1'].sum()
            ra['c2'] = ra['c2'] / ra['c2'].sum()
        else:
            ra['c1'] = ra['c1'] * randnrs3
            ra['c2'] = ra['c2'] * randnrs3
            ra['c1'] = ra['c1'] / ra['c1'].sum()
            ra['c2'] = ra['c2'] / ra['c2'].sum()

        df['c1'] = ra['c1'] * samplecount
        df['c2'] = ra['c2'] * samplecount
        df = df.applymap(round)
        for c in ['c1', 'c2']:
            diff = int(samplecount - df[c].sum())
            if diff > 0:
                indx = list(df[c].sample(n=diff, weights=ra[c]).index)
                df.loc[indx, c] = df.loc[indx, c] + 1
            elif diff < 0:
                indx = list(df[df[c] > 0].index)
                np.random.shuffle(indx)
                df.loc[indx[:-diff], c] = df.loc[indx[:-diff], c]-1

        dtemp = df.copy()
        dtemp['prevc1'] = prevc1
        dtemp['prevc2'] = prevc2

        alfa0 = naiveDivAlpha(dtemp, 0)
        alfa1 = naiveDivAlpha(dtemp, 1)
        alfa2 = naiveDivAlpha(dtemp, 2)
        beta0 = naiveDivBeta(dtemp, 0)
        beta1 = naiveDivBeta(dtemp, 1)
        beta2 = naiveDivBeta(dtemp, 2)
        dist0 = beta2Dist(beta0, 0, 'naive')
        dist1 = beta2Dist(beta1, 1, 'naive')
        dist2 = beta2Dist(beta2, 2, 'naive')

        outi = i
        dt0 = np.average([dist0.loc['c1', 'prevc1'], dist0.loc['c2', 'prevc2']])
        dt1 = np.average([dist1.loc['c1', 'prevc1'], dist1.loc['c2', 'prevc2']])
        dt2 = np.average([dist2.loc['c1', 'prevc1'], dist2.loc['c2', 'prevc2']])
        dx0 = dist0.loc['c1', 'c2']
        dx1 = dist1.loc['c1', 'c2']
        dx2 = dist2.loc['c1', 'c2']
        a0 = np.average(alfa0.values)
        a1 = np.average(alfa1.values)
        a2 = np.average(alfa2.values)
        output.append([outi, dt0, dt1, dt2, dx0, dx1, dx2, a0, a1, a2])
        prevc1 = df['c1']
        prevc2 = df['c2']

        if i % 50 == 0:
            print(i, '  ', diff)

    outdf = pd.DataFrame(output[1:], columns=output[0])
    outdf.to_csv('Deterministicmodel.csv')

    plt.figure(figsize=(10, 8))
    plt.plot(outdf['iter'], outdf['dt_0ND'], label='dt_0ND')
    plt.plot(outdf['iter'], outdf['dt_1ND'], label='dt_1ND')
    plt.plot(outdf['iter'], outdf['dt_2ND'], label='dt_2ND')
    plt.legend(loc='best')
    plt.show()

    plt.figure(figsize=(10, 8))
    plt.plot(outdf['iter'], outdf['dx_0ND'], label='dx_0ND')
    plt.plot(outdf['iter'], outdf['dx_1ND'], label='dx_1ND')
    plt.plot(outdf['iter'], outdf['dx_2ND'], label='dx_2ND')
    plt.legend(loc='best')
    plt.show()

    plt.figure(figsize=(10, 8))
    plt.plot(outdf['iter'], outdf['a_0ND'], label='a_0ND')
    plt.plot(outdf['iter'], outdf['a_1ND'], label='a_1ND')
    plt.plot(outdf['iter'], outdf['a_2ND'], label='a_2ND')
    plt.legend(loc='best')
    plt.show()


# Raup-Crick according to Stegen 2013 (https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r)
# Compared null simulation to beta diversities based on Hill (q=0,1,2) and Bray-Curtis
# Returns dictionary with data matrices as Hill_0, Hill_1, Hill_2, and Bray
def raup(tab, iters=999):
    Hillreal0 = naiveDivBeta(tab, 0)
    Hillreal1 = naiveDivBeta(tab, 1)
    Hillreal2 = naiveDivBeta(tab, 2)
    BCreal = bray(tab)

    output0 = pd.DataFrame(0, index=tab.columns, columns=tab.columns)
    output1 = pd.DataFrame(0, index=tab.columns, columns=tab.columns)
    output2 = pd.DataFrame(0, index=tab.columns, columns=tab.columns)
    outputBC = pd.DataFrame(0, index=tab.columns, columns=tab.columns)

    dfBin = tab.copy()
    dfBin[dfBin > 0] = 1
    frequency = dfBin.sum(axis=1)  # Gets SV frequency across samples
    abundance = tab.sum(axis=1)  # Get total count frequency

    dfOnes = pd.DataFrame(1, index=tab.index, columns=tab.columns)
    for i in range(iters):
        print('iteration= ', i, ' out of ', iters)
        dfNull = pd.DataFrame(0, index=tab.index, columns=tab.columns)

        for smp in tab.columns:
            # Subsample SVs based on frequency
            subsampleSVs = list(dfOnes[smp].sample(dfBin[smp].sum(), replace=False, weights=frequency).index)
            dfNull.loc[subsampleSVs, smp] = 1
            # Assign counts to SVs based on abundance
            subsampleCounts = dfOnes.loc[subsampleSVs, smp].sample(tab[smp].sum() - dfNull[smp].sum(), replace=True,
                                                                   weights=abundance)
            subsampleCounts = subsampleCounts.groupby(subsampleCounts.index).sum()
            dfNull[smp] = dfNull[smp].add(subsampleCounts, fill_value=0)
        # Calculate Beta diversities
        sim0 = naiveDivBeta(dfNull, 0)
        sim1 = naiveDivBeta(dfNull, 1)
        sim2 = naiveDivBeta(dfNull, 2)
        simBC = bray(dfNull)

        # Add 1 to output matrix if observed dissimilarities are larger than simulated
        output0[Hillreal0 > sim0] = output0[Hillreal0 > sim0] + 1
        output0[Hillreal0 == sim0] = output0[Hillreal0 == sim0] + 0.5
        output1[Hillreal1 > sim1] = output1[Hillreal1 > sim1] + 1
        output1[Hillreal1 == sim1] = output1[Hillreal1 == sim1] + 0.5
        output2[Hillreal2 > sim2] = output2[Hillreal2 > sim2] + 1
        output2[Hillreal2 == sim2] = output2[Hillreal2 == sim2] + 0.5
        outputBC[BCreal > simBC] = outputBC[BCreal > simBC] + 1
        outputBC[BCreal == simBC] = outputBC[BCreal == simBC] + 0.5

    # Divide by number of iterations to get fractions in range -1 to 1
    output0 = output0 / iters
    output0 = (output0-0.5)*2
    output1 = output1/iters
    output1 = (output1-0.5)*2
    output2 = output2/iters
    output2 = (output2-0.5)*2
    outputBC = outputBC/iters
    outputBC = (outputBC-0.5)*2
    out = {}
    out['Hill_0'] = output0
    out['Hill_1'] = output1
    out['Hill_2'] = output2
    out['Bray'] = outputBC
    return out


def bNTI(tab, distmat, iters=999):
    ntdreal = bMNTD(tab, distmat)

    output = np.zeros((len(tab.columns), len(tab.columns), iters))
    for i in range(iters):
        print('Iteration ', i, ' out of ', iters)
        svs = list(distmat.index);
        np.random.shuffle(svs)
        dmrand = pd.DataFrame(distmat.values, index=svs, columns=svs)
        ntdrand = bMNTD(tab, dmrand)
        output[:, :, i] = ntdrand

    outputmean = np.mean(output, axis=2)
    outputstdev = np.std(output, axis=2)
    ses = (ntdreal - outputmean) / outputstdev
    return ses
