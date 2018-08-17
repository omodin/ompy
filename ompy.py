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
import prince
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import seaborn as sns
from scipy import optimize


pd.options.mode.chained_assignment = None  # default='warn'

# Default colors used in plotting
colorlist1 = ['black', 'cyan', 'blue', 'blue', 'blue', 'red', 'darkred', 'darkred', 'darkred', 'grey',
             'darkgrey', 'darkgrey', 'darkgrey', 'gold', 'orange', 'orange', 'orange', 'darkmagenta', 'magenta',
              'fuchsia', 'pink', 'Burlywood', 'green', 'chartreuse', 'Burlywood']

colorlistX = ['black', 'darkblue', 'darkred', 'gold', 'darkgrey', 'blue', 'red', 'orange', 'grey', 'darkblue', 'brown',
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
    if path == 'None':
        print('Error: No path specified')
        return 0

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

# Function that makes sure different objects have the same SV names
def calibrateSVs(objlist):
    #Give all SVs new names, associate sequence with name
    svlist = []
    for i in range(len(objlist)):
        svlist = svlist + objlist[i]['seq']['seq'].tolist()
    svlist = list(set(svlist))
    svdict = {}; counter = 0
    for sv in svlist:
        counter += 1
        svdict[sv] = 'SV'+str(counter)

    #Change the name of all SVs
    for i in range(len(objlist)):
        seq = objlist[i]['seq']
        seq['newSV'] = np.nan
        if 'tab' in objlist[i].keys():
            tab = objlist[i]['tab']
            ra = objlist[i]['ra']
            tab['newSV'] = np.nan
            ra['newSV'] = np.nan
        if 'tax' in objlist[i].keys():
            tax = objlist[i]['tax']
            tax['newSV'] = np.nan

        for n in seq.index:
            newSVname = svdict[seq.loc[n, 'seq']]
            seq.loc[n, 'newSV'] = newSVname
            if 'tab' in objlist[i].keys():
                tab.loc[n, 'newSV'] = newSVname
                ra.loc[n, 'newSV'] = newSVname
            if 'tax' in objlist[i].keys():
                tax.loc[n, 'newSV'] = newSVname

        seq = seq.set_index('newSV')
        if 'tab' in objlist[i].keys():
            tab = tab.set_index('newSV')
            ra = ra.set_index('newSV')
        if 'tax' in objlist[i].keys():
            tax = tax.set_index('newSV')

        objlist[i]['seq'] = seq
        if 'tab' in objlist[i].keys():
            objlist[i]['tab'] = tab
            objlist[i]['ra'] = ra
        if 'tax' in objlist[i].keys():
            objlist[i]['tax'] = tax
    return objlist

# Function that merges different objects
def combineObjects(objlist):
    calibratedObjects = calibrateSVs(objlist)

    #Change sample names to include an index indicating object
    for i in range(len(calibratedObjects)):
        colnames = calibratedObjects[i]['tab'].columns.tolist()
        newcolnames = []
        for n in colnames:
            newcolnames.append('Obj'+str(i)+'_'+n)
        calibratedObjects[i]['tab'].columns = newcolnames
        calibratedObjects[i]['ra'].columns = newcolnames

        rownames = calibratedObjects[i]['meta'].index.tolist()
        newrownames = []
        for n in rownames:
            newrownames.append('Obj'+str(i)+'_'+n)
        calibratedObjects[i]['meta']['newindex'] = newrownames
        calibratedObjects[i]['meta'] = calibratedObjects[i]['meta'].set_index('newindex')

    #Merge all dataframes
    joined_object = calibratedObjects[0].copy()
    for i in range(1, len(calibratedObjects)):
        tab0 = joined_object['tab']
        tab1 = calibratedObjects[i]['tab']
        joined_object['tab'] = tab0.join(tab1, how='outer').fillna(0)

        ra0 = joined_object['ra']
        ra1 = calibratedObjects[i]['ra']
        joined_object['ra'] = ra0.join(ra1, how='outer').fillna(0)

        SV0list = joined_object['seq'].index.tolist()
        SV1list = calibratedObjects[i]['seq'].index.tolist()
        addto0 = [x for x in SV1list if x not in SV0list]
        joined_object['seq'] = pd.concat([joined_object['seq'], calibratedObjects[i]['seq'].loc[addto0, :]], join='outer')

        if 'tax' in calibratedObjects[i].keys():
            tx0list = joined_object['tax'].index.tolist()
            tx1list = calibratedObjects[i]['tax'].index.tolist()
            addto0 = [x for x in tx1list if x not in tx0list]
            joined_object['tax'] = pd.concat([joined_object['tax'], calibratedObjects[i]['tax'].loc[addto0, :]], join='outer')

        meta0 = joined_object['meta']
        meta1 = calibratedObjects[i]['meta']
        joined_object['meta'] = pd.concat([meta0, meta1], join='outer')

    joined_object['tab'] = joined_object['tab'].applymap(int)
    return joined_object


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

# Funtion that subset dataframes based on a certain number of the most abundant SVs
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

# Subset obj based on text in taxa
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

# Merges samples based on variable and list from meta data
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

# Table with rarefied counts based on calculation, similar to rarafaction without replacement
def rarefy1(tab, depth='min'):
    tab = tab.applymap(int) #Make sure table elements are integers
    if depth == 'min':  # Define read depth
        reads = min(tab.sum())
    else:
        reads = depth

    samples = tab.columns.tolist()
    svs = tab.index.tolist()
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
            addSV = tab[smp].sample(n=addNR, replace=False, weights=diffs).index.tolist()
            avrundad[addSV] = avrundad[addSV] + 1
        rtab[smp] = avrundad
    return rtab

# Table with rarefied counts, Randomly pick depth counts of present SVs.
# Chance of being picked is proportional to number of counts (i.e. with replacement).
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

# Plots heatmap, order is the heading in metadata for the column that specifies logical order of samples
def plotHeatmap(obj, var='None', levels=['Phylum', 'Genus'], subsetLevels='None', subsetPatterns='None',
                order='None', numberToPlot=20, method='max_sample', nameType='SV', labels=True, labelSize=1, fontSize=15, saveName='None'):
    #Merge samples based on var
    if var != 'None':
        merged_obj = mergeSamples(obj, var=var)
    else:
        merged_obj = obj.copy()

    #Calculate relative abundances and store in df ra
    tab = merged_obj['tab']
    ra = 100*tab/tab.sum()

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
        ra = subset_obj['ra']

    ## Groupby taxa
    obj_to_group = {}
    obj_to_group['ra'] = ra; obj_to_group['tax'] = merged_obj['tax']
    taxa_obj = groupbyTaxa(obj_to_group, levels=levels, nameType=nameType)
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

    elif method == 'max_sample_all':
        retain = []
        for c in cols:
            ra = ra.sort_values(by=[c], ascending=False)
            svs = ra.index.tolist()[:number]
            retain = retain + svs
        retain = list(set(retain))

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
                    labelvalues.loc[r, c] = -0.1
                elif value < 10 and value >= 0.1:
                    labelvalues.loc[r, c] = round(value, 1)
                elif value > 99:
                    labelvalues.loc[r, c] = int(99)
                elif value >= 10:
                    labelvalues.loc[r, c] = round(value)
                else:
                    labelvalues.loc[r, c] = 0

    #Plot
    fig, ax = plt.subplots(figsize=(14, 10))
    sns.set(font_scale=labelSize)
    if labels:
        sns.heatmap(table, annot=labelvalues, cmap='Reds', linewidths=0.5, robust=True, cbar=False, ax=ax)
    else:
        sns.heatmap(table, cmap='Reds', linewidths=0.5, robust=True, cbar=False, ax=ax)
    plt.xticks(rotation=90)
    plt.setp(ax.get_xticklabels(), fontsize=fontSize)
    ax.set_ylabel('')
    plt.setp(ax.get_yticklabels(), fontsize=fontSize)
    plt.tight_layout()
    if saveName != 'None':
        plt.savefig(saveName)
    plt.show()

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
def plotPCoA(dist, meta, var1='None', var2='None', tag='None', order='None', title='', colorlist=colorlistX, markerlist=markerlist, savename='None'):
    # Do the PCoA in ecopy and get results for PC1 and PC2 into dataframe
    pcs = ecopy.pcoa(dist, correction='1')
    pc_values = ecopy.pcoa.summary(pcs)
    pc1_perc = round(100 * pc_values.iloc[1, 0], 1)
    xn = 'PC1 (' + str(pc1_perc) + '%)'
    pc2_perc = round(100 * pc_values.iloc[1, 1], 1)
    yn = 'PC2 (' + str(pc2_perc) + '%)'

    koord = pcs.biplot(coords=True)['Objects']
    smplist = dist.index
    pcoadf = pd.DataFrame(koord, index=smplist, columns=[xn, yn])

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

    # Do the plotting
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
                    ax.scatter(xlist, ylist, label=None, color=colorlist[i], marker=markerlist[jcounter], s=70)

                    if jcounter not in shapeTracker:
                        linesShape[0].append(ax.scatter([], [], label=str(smpcats2[j]), color='black', marker=markerlist[jcounter]))
                        linesShape[1].append(smpcats2[jcounter])
                        shapeTracker.append(jcounter)
                jcounter += 1

            # Here set both legends for color and marker
            ax.legend(linesColor[0], linesColor[1], ncol=1, bbox_to_anchor=(1.3, 1), title='Time', frameon=False, markerscale=1.5, fontsize=18)
            from matplotlib.legend import Legend
            leg = Legend(ax, linesShape[0], linesShape[1], bbox_to_anchor=(1.3, 0.4), title='Reactor', frameon=False, markerscale=1.5, fontsize=18)
            ax.add_artist(leg)

        else: #If there is no var2, change both color and marker with each category in var1
            linesColor[0].append(ax.scatter([], [], label=str(smpcats1[i]), color=colorlist[i], marker=markerlist[i]))
            linesColor[1].append(smpcats1[i])

            xlist = metaPlot_i[xn]
            ylist = metaPlot_i[yn]
            ax.scatter(xlist, ylist, label=None, color=colorlist[i], marker=markerlist[i], s=100)
            ax.legend(linesColor[0], linesColor[1], ncol=1, bbox_to_anchor=(1.0, 1.05), title='', frameon=False, markerscale=1.5, fontsize=12)
    # for lgi in range(len(smpcats1)):
    #     lgnd.legendHandles[lgi]._sizes = [120]

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

def plotMFA(objlist, nameslist, path, rarefy=False):
    if len(objlist) != len(nameslist):
        print('Lengths dont match')
        return 0

    #Merge tabs into new df, make sure SVs/OTUs are common
    df = objlist[0]['tab']
    newcolnames = [] #Set new column names to indicate the object
    [newcolnames.append(item + ':' + nameslist[0]) for item in df.columns]
    df.columns = newcolnames
    seq = objlist[0]['seq']
    df = pd.concat([df, seq], axis=1)
    df = df.set_index('seq')
    for i in range(1, len(objlist)):
        tab = objlist[i]['tab']
        newcolnames = []  # Set new column names to indicate the object
        [newcolnames.append(item + ':' + nameslist[i]) for item in tab.columns]
        tab.columns = newcolnames
        seq = objlist[i]['seq']
        tab = pd.concat([tab, seq], axis=1)
        tab = tab.set_index('seq')
        df = pd.concat([df, tab], axis=1)
    df = df.fillna(0)

    seqlist = []; seqnames = []; counter = 0
    for s in df.index:
        seqlist.append(s)
        counter += 1
        seqnames.append('SV'+str(counter))
    df['SV'] = seqnames
    df = df.set_index('SV') #New df with all tabs

    if rarefy:
        df = rarefy1(df)

    #Prepare df for MFA - prepare group dictionary and df
    smp_dict = {}
    counter = 0
    for c in df.columns:
        smp = c.split(':')[0]
        if smp not in smp_dict:
            list(range(counter * len(df.index), counter * len(df.index) + len(df.index)))
            smp_dict[smp] = [str(x) for x in range(counter*len(df.index), counter*len(df.index) + len(df.index))]
            counter += 1
    mfa_df = pd.DataFrame(0, index=nameslist, columns=[str(x) for x in range(len(df.index)*len(smp_dict))])

    #Prepare df for MFA - put values in df
    for c in df.columns:
        print(c)
        smp = c.split(':')[0]
        tb = c.split(':')[1]
        mfa_df.loc[tb, smp_dict[smp]] = df[c].values.tolist()

    #Remove columns with all zeros
    keepcols = mfa_df.sum()[mfa_df.std() != 0].index.tolist()
    remcols = mfa_df.sum()[mfa_df.std() == 0].index.tolist()
    for smp in smp_dict.keys():
        print(smp)
        dict_keepscols = smp_dict[smp]
        intersect = list(set(remcols).intersection(dict_keepscols))
        for item in intersect:
            remcols.remove(item)
            dict_keepscols.remove(item)
        smp_dict[smp] = dict_keepscols
    mfa_df = mfa_df[keepcols]

    #Save data
    if rarefy:
        mfa_df.to_csv(path+'MFA_df_rarefied.csv')
        with open(path+'MFA_groups_rarefied.pickle', 'wb') as f:
            pickle.dump(smp_dict, f)
    else:
        mfa_df.to_csv(path+'MFA_df.csv')
        with open(path+'MFA_groups.pickle', 'wb') as f:
            pickle.dump(smp_dict, f)
    print('Data saved to path')

    ## Do MFA
    mfa = prince.MFA(groups=smp_dict, n_components=2, n_iter=3, copy=True, engine='auto', random_state=42)
    mfa = mfa.fit(mfa_df)
    print(mfa.row_coordinates(mfa_df))
    mfa.plot_row_coordinates(mfa_df)
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

# Compare different empirical functions as fits to Diversity Accumulation Curves
# Returns a dataframe wih r2 and parameters for each function and curve
# If savename is specified, a table of the rarefaction curves and each plot are saved
def compareDAC(obj, var, slist='all', step=2000, q=0, savename='None'):
    tab = obj['tab']
    meta = obj['meta']

    if slist == 'all':
        slist = meta[var]

    # Pick out the samples in tab that are specified in slist
    submeta = meta.loc[meta[var].isin(slist)]
    if 'Logic_order' in submeta.columns:
        submeta = submeta.sort_values(by=['Logic_order']) #Make sure samples  in smplist are in correct order
    smplist = submeta.index.tolist()
    subtab = tab.loc[:, smplist]
    subtab['sum'] = subtab.sum(axis=1)
    subtab = subtab[subtab['sum'] > 0]
    subtab = subtab.drop('sum', axis=1)

    plottab = [['count']+smplist] #Table to plot, first col is count
    plottab.append([1]*(len(smplist)+1))
    maxcount = int(max(subtab.sum())) #Max number of reads in a sample
    numberofsteps = math.floor(maxcount/step) #Number of steps to take
    counter = 0  #Just to keep track of progress
    # Start iteration
    for i in range(step, maxcount, step):
        counter+=1
        if counter%5 == 0:  #Just to keep track of progress
            print(counter, ' out of ', numberofsteps + 1, ' steps')

        # Rarefy to specific depth
        raretab = rarefy1(subtab, depth=i)
        svcount = naiveDivAlpha(raretab, q)  #Calculate diversity value (order q)
        templist=[i]
        for s in smplist:
            if s in svcount.index:
                templist = templist + [svcount[s]]
            else:
                templist = templist + [0]
        plottab.append(templist)
    plottab = pd.DataFrame(plottab[1:], columns=plottab[0])

    # Data smoothing
    # Flather 1996 J Biogeography, 23, 155; Tjörve 2003 J Biogeography, 30, 827;
    def monod(xdata, a, b):
        return a*xdata/(b+xdata)
    monod_dict = {}
    def power(xdata, a, b):
        return a*xdata**b
    power_dict = {}
    def exponential(xdata, a, b):
        return a+b*np.log(xdata)
    exponential_dict = {}
    def neg_exponential(xdata, a, b):
        return a*(1-np.exp(-b*xdata))
    neg_exponential_dict = {}
    def asym_regression(xdata, a, b, c):
        return a-b*c**xdata
    asym_regression_dict = {}
    def rational(xdata, a, b, c):
        return (a+b*xdata)/(1+c*xdata)
    rational_dict = {}

    # Help function
    def closest_value_pos(arr, val): #Find position of a value in array that is closest to val
        return np.abs(arr-val).argmin()
    def r_squared(ydata, fdata): # Find r squared value for two lists of data of equal length
        y_mean = np.mean(ydata)
        ss_tot = 0
        ss_res = 0
        for d in range(len(ydata)):
            ss_tot = ss_tot + (ydata[d] - y_mean)**2
            ss_res = ss_res + (ydata[d] - fdata[d])**2
        return 1 - ss_res/ss_tot

    # Use plottab2 for data smoothing
    plottab = plottab.applymap(float)
    for smp in smplist: #Start from 1 because 0 is 'count'
        xydata = plottab.loc[:, ['count', smp]][plottab.loc[:, smp] != 0]

        #Monod
        a_est = xydata.iloc[:, 1].tolist()[-1] #Final value of y is estimate of a_est
        b_pos = closest_value_pos(xydata.iloc[:, 1], a_est/2) #Position of half saturation constant
        b_est = xydata.iloc[:, 0].tolist()[b_pos] #Estimate of ks
        popt, pcov =optimize.curve_fit(monod, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est])
        monod_dict[smp] = popt #These are fitted Monod constants

        #Power
        a_est = 100
        b_est = 0.1
        popt, pcov = optimize.curve_fit(power, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est])
        power_dict[smp] = popt

        #Exponential
        a_est = -100
        b_est = 100
        popt, pcov = optimize.curve_fit(exponential, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est])
        exponential_dict[smp] = popt

        #Negative exponential
        a_est = 1000
        b_est = 10**-8
        popt, pcov = optimize.curve_fit(neg_exponential, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est])
        neg_exponential_dict[smp] = popt

        #Asymptotic regression
        a_est = 1000
        b_est = 500
        c_est = 0.1
        popt, pcov = optimize.curve_fit(asym_regression, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est, c_est])
        asym_regression_dict[smp] = popt

        #Rational
        a_est = 500
        b_est = 1
        c_est = 0.001
        popt, pcov = optimize.curve_fit(rational, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est, c_est])
        rational_dict[smp] = popt

        ###perr = np.sqrt(np.diag(pcov))  #This can be used to find standard deviations of constants

    # Save results as plottab and print as csv file if desired
    if savename != 'None':
        plottab.to_csv(savename+'.csv')

    #Run optimization and plot results
    regression_output = [['Function', 'Sample', 'r_squared', 'Constants']]
    plt.rcParams.update({'font.size': 16})
    for func in [[monod, monod_dict, 'Monod'], [power, power_dict, 'Power'], [exponential, exponential_dict, 'Exp'],
                 [neg_exponential, neg_exponential_dict, 'Neg_exp'], [asym_regression, asym_regression_dict, 'Asym_regr'],
                 [rational, rational_dict, 'Rational']]:

        colorcheck = []
        plt.figure(figsize=(10, 6))

        for smp in smplist:
            stype = submeta.loc[smp, var]
            xdata = plottab['count'][plottab[smp] != 0]
            ydata = plottab[smp][plottab[smp] != 0]
            fdata = func[0](xdata, *func[1][smp])

            regression_data = [func[2], smp, r_squared(ydata, fdata),  func[1][smp]. tolist()]
            regression_output.append(regression_data)

            #Plot
            if stype not in colorcheck:
                colorcheck.append(stype)
                plt.plot(xdata, ydata, lw=0, marker='.', ms=2, label='_nolegend_', color=colorlistX[len(colorcheck)-1])
                plt.plot(xdata, func[0](xdata, *func[1][smp]), lw=1, label=stype, color=colorlistX[len(colorcheck)-1])
            else:
                plt.plot(xdata, ydata, lw=0, marker='.', ms=2, label='_nolegend_', color=colorlistX[len(colorcheck) - 1])
                plt.plot(xdata, func[0](xdata, *func[1][smp]), lw=1, label='_nolegend_', color=colorlistX[len(colorcheck) - 1])
        plt.legend(loc='best')
        plt.xlabel('Number of reads')
        plt.xticks(rotation=90)
        plt.ylabel('Diversity')
        plt.title(func[2])
        plt.tight_layout()
        if savename != 'None':
            plt.savefig(savename+'_'+func[3]+'.png')
        plt.show()

    regression_output = pd.DataFrame(regression_output[1:], columns=regression_output[0])
    return regression_output

# Plots smoothed Diversity Accumulation Curves based on chosen function
# Also prints R2 value for each line
def plotSmoothedDAC(obj, var, slist='all', step=2000, q=0, function='rational', savename='None'):
    tab = obj['tab']
    meta = obj['meta']

    if slist == 'all':
        slist = meta[var]

    # Pick out the samples in tab that are specified in slist
    submeta = meta.loc[meta[var].isin(slist)]
    if 'Logic_order' in submeta.columns:
        submeta = submeta.sort_values(by=['Logic_order']) #Make sure samples  in smplist are in correct order
    smplist = submeta.index.tolist()
    subtab = tab.loc[:, smplist]
    subtab['sum'] = subtab.sum(axis=1)
    subtab = subtab[subtab['sum'] > 0]
    subtab = subtab.drop('sum', axis=1)

    plottab = [['count']+smplist] #Table to plot, first col is count
    plottab.append([1]*(len(smplist)+1))
    maxcount = int(max(subtab.sum())) #Max number of reads in a sample
    numberofsteps = math.floor(maxcount/step) #Number of steps to take
    counter = 0  #Just to keep track of progress
    # Start iteration
    for i in range(step, maxcount, step):
        counter+=1
        if counter%5 == 0:  #Just to keep track of progress
            print(counter, ' out of ', numberofsteps + 1, ' steps')

        # Rarefy to specific depth
        raretab = rarefy1(subtab, depth=i)
        svcount = naiveDivAlpha(raretab, q)  #Calculate diversity value (order q)
        templist=[i]
        for s in smplist:
            if s in svcount.index:
                templist = templist + [svcount[s]]
            else:
                templist = templist + [0]
        plottab.append(templist)
    plottab = pd.DataFrame(plottab[1:], columns=plottab[0])

    # Data smoothing
    # Flather 1996 J Biogeography, 23, 155; Tjörve 2003 J Biogeography, 30, 827;
    def monod(xdata, a, b):
        return a*xdata/(b+xdata)
    def power(xdata, a, b):
        return a*xdata**b
    def exponential(xdata, a, b):
        return a+b*np.log(xdata)
    def neg_exponential(xdata, a, b):
        return a*(1-np.exp(-b*xdata))
    def asym_regression(xdata, a, b, c):
        return a-b*c**xdata
    def rational(xdata, a, b, c):
        return (a+b*xdata)/(1+c*xdata)

    # Help function
    def closest_value_pos(arr, val): #Find position of a value in array that is closest to val
        return np.abs(arr-val).argmin()
    def r_squared(ydata, fdata): # Find r squared value for two lists of data of equal length
        y_mean = np.mean(ydata)
        ss_tot = 0
        ss_res = 0
        for d in range(len(ydata)):
            ss_tot = ss_tot + (ydata[d] - y_mean)**2
            ss_res = ss_res + (ydata[d] - fdata[d])**2
        return 1 - ss_res/ss_tot

    # Use plottab2 for original data, put smoothed data back into plottab
    plottab2 = plottab.copy()
    plottab2 = plottab2.applymap(float)
    par_dict = {}
    for smp in smplist:
        xydata = plottab2.loc[:, ['count', smp]][plottab2.loc[:, smp] != 0]

        if function == 'rational':
            a_est = 500
            b_est = 1
            c_est = 0.001
            popt, pcov = optimize.curve_fit(rational, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est, c_est])
            plottab.loc[xydata.index, smp] = rational(xydata.iloc[:, 0], *popt)
            par_dict[smp] = popt; plotfunc = rational

        elif function == 'monod':
            a_est = xydata.iloc[:, 1].tolist()[-1] #Final value of y is estimate of a_est
            b_pos = closest_value_pos(xydata.iloc[:, 1], a_est/2) #Position of half saturation constant
            b_est = xydata.iloc[:, 0].tolist()[b_pos] #Estimate of ks
            popt, pcov =optimize.curve_fit(monod, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est])
            plottab.loc[xydata.index, smp] = monod(xydata.iloc[:, 0], *popt)
            par_dict[smp] = popt; plotfunc = monod

        elif function == 'power':
            a_est = 100
            b_est = 0.1
            popt, pcov = optimize.curve_fit(power, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est])
            plottab.loc[xydata.index, smp] = power(xydata.iloc[:, 0], *popt)
            par_dict[smp] = popt; plotfunc = power

        elif function == 'exponential':
            a_est = -100
            b_est = 100
            popt, pcov = optimize.curve_fit(exponential, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est])
            plottab.loc[xydata.index, smp] = exponential(xydata.iloc[:, 0], *popt)
            par_dict[smp] = popt; plotfunc = exponential

        elif function == 'neg_exponential':
            a_est = 1000
            b_est = 10**-8
            popt, pcov = optimize.curve_fit(neg_exponential, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est])
            plottab.loc[xydata.index, smp] = neg_expontential(xydata.iloc[:, 0], *popt)
            par_dict[smp] = popt; plotfunc = neg_exponential

        elif function == 'asym_regression':
            a_est = 1000
            b_est = 500
            c_est = 0.1
            popt, pcov = optimize.curve_fit(asym_regression, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est, c_est])
            plottab.loc[xydata.index, smp] = asym_regression(xydata.iloc[:, 0], *popt)
            par_dict[smp] = popt; plotfunc = asym_regression

        else:
            return 'Error, function not found'

    # Print smoothed results as csv file if desired
    if savename != 'None':
        plottab.to_csv(savename+'.csv')

    #Plot smoothed results
    plt.rcParams.update({'font.size': 16})
    colorcheck = []
    plt.figure(figsize=(10, 6))

    for smp in smplist:
        fdata = plottab[smp][plottab[smp] != 0]
        ydata = plottab2[smp][plottab2[smp] != 0]
        print(smp,' _R2= ', r_squared(ydata, fdata))
        print(*par_dict[smp])

        totcount = subtab[smp].sum()
        xdata2 = range(10, totcount, 10)
        fdata2 = plotfunc(xdata2, *par_dict[smp])

        stype = submeta.loc[smp, var]
        if stype not in colorcheck:
            colorcheck.append(stype)
            plt.plot(xdata2, fdata2, lw=1, label=stype, color=colorlistX[len(colorcheck)-1])
        else:
            plt.plot(xdata2, fdata2, lw=1, label='_nolegend_', color=colorlistX[len(colorcheck) - 1])
    plt.legend(loc='best')
    plt.xlabel('Number of reads')
    plt.xticks(rotation=90)
    plt.ylabel('Diversity')
    plt.tight_layout()
    if savename != 'None':
        plt.savefig(savename+'_'+func[3]+'.png')
    plt.show()

def plotMonodRationalDAC(obj, var, slist='all', step=2000, q=0, savename='None', plot='Yes',
                         colorlist='None', xlim='None', ylim='None', fontsize=22, figsize=(10, 6), printprogress='Yes'):
    tab = obj['tab']
    meta = obj['meta']

    if slist == 'all':
        slist = meta[var]

    # Pick out the samples in tab that are specified in slist
    submeta = meta.loc[meta[var].isin(slist)]
    if 'Logic_order' in submeta.columns:
        submeta = submeta.sort_values(by=['Logic_order']) #Make sure samples  in smplist are in correct order
    smplist = submeta.index.tolist()
    subtab = tab.loc[:, smplist]
    subtab['sum'] = subtab.sum(axis=1)
    subtab = subtab[subtab['sum'] > 0]
    subtab = subtab.drop('sum', axis=1)

    #Prepare for initial iteration
    plottab_low = [['count']+smplist] #Table with only low counts, first col is count
    mincount = int(min(subtab.sum())) #Min number of reads in a sample
    if mincount > step:
        mincount = step
    shortsteps = int(mincount/50)
    listofcounts_low = [10]+list(range(shortsteps, mincount, shortsteps))

    #Prepare for full iteration
    plottab = [['count']+smplist] #Full table to plot, first col is count
    maxcount = int(max(subtab.sum())) #Max number of reads in a sample
    maxlist = list(set(subtab.sum().values.tolist()))
    rangelist = list(range(step, maxcount, step))
    listofcounts = maxlist + rangelist
    listofcounts = sorted(listofcounts)

    counter = 0  #Just to keep track of progress

    # Start for low count iteration
    for i in listofcounts_low:
        if printprogress == 'Yes':
            counter+=1
            if counter%10 == 0:  #Just to keep track of progress
                print(counter, ' out of ', len(listofcounts_low)+len(listofcounts), ' steps')

        # Rarefy to specific depth
        raretab = rarefy1(subtab, depth=i)
        svcount = naiveDivAlpha(raretab, q)  #Calculate diversity value (order q)
        templist=[i]
        for s in smplist:
            if s in svcount.index:
                templist = templist + [svcount[s]]
            else:
                templist = templist + [0]
        plottab_low.append(templist)
    plottab_low = pd.DataFrame(plottab_low[1:], columns=plottab_low[0])

    # Start for full count iteration
    for i in listofcounts:
        if printprogress == 'Yes':
            counter+=1
            if counter%10 == 0:  #Just to keep track of progress
                print(counter, ' out of ', len(listofcounts_low) + len(listofcounts), ' steps')

        # Rarefy to specific depth
        raretab = rarefy1(subtab, depth=i)
        svcount = naiveDivAlpha(raretab, q)  #Calculate diversity value (order q)
        templist=[i]
        for s in smplist:
            if s in svcount.index:
                templist = templist + [svcount[s]]
            else:
                templist = templist + [0]
        plottab.append(templist)
    plottab = pd.DataFrame(plottab[1:], columns=plottab[0])

    # Data smoothing
    # Flather 1996 J Biogeography, 23, 155; Tjörve 2003 J Biogeography, 30, 827;
    def monod(xdata, a, b):
        return a*xdata/(b+xdata)
    parameter_dict_low = {}
    def rational(xdata, a, b, c):
        return (a+b*xdata)/(1+c*xdata)
    parameter_dict = {}

    # Help function
    def closest_value_pos(arr, val): #Find position of a value in array that is closest to val
        return np.abs(arr-val).argmin()
    def r_squared(ydata, fdata): # Find r squared value for two lists of data of equal length
        y_mean = np.mean(ydata)
        ss_tot = 0
        ss_res = 0
        for d in range(len(ydata)):
            ss_tot = ss_tot + (ydata[d] - y_mean)**2
            ss_res = ss_res + (ydata[d] - fdata[d])**2
        return 1 - ss_res/ss_tot

    # Put smoothed data back into plottab_low and plottab
    orig_plottab_low = plottab_low.copy()
    orig_plottab = plottab.copy()

    plottab_low = plottab_low.applymap(float)
    plottab = plottab.applymap(float)
    for smp in smplist:
        #Monod for low
        xydata_low = plottab_low.loc[:, ['count', smp]][plottab_low.loc[:, smp] != 0]
        a_est = xydata_low.iloc[:, 1].tolist()[-1] #Final value of y is estimate of a_est
        b_pos = closest_value_pos(xydata_low.iloc[:, 1], a_est/2) #Position of half saturation constant
        b_est = xydata_low.iloc[:, 0].tolist()[b_pos] #Estimate of ks
        popt, pcov =optimize.curve_fit(monod, xydata_low.iloc[:, 0], xydata_low.iloc[:, 1], p0=[a_est, b_est])
        plottab_low.loc[xydata_low.index, smp] = monod(xydata_low.iloc[:, 0], *popt)
        parameter_dict_low[smp] = popt

        #Rational for high
        xydata = plottab.loc[:, ['count', smp]][plottab.loc[:, smp] != 0]
        a_est = 500
        b_est = 1
        c_est = 0.001
        popt, pcov = optimize.curve_fit(rational, xydata.iloc[:, 0], xydata.iloc[:, 1], p0=[a_est, b_est, c_est])
        plottab.loc[xydata.index, smp] = rational(xydata.iloc[:, 0], *popt)
        parameter_dict[smp] = popt

    #Combine plottab_low and plottab
    plottab = plottab[plottab['count'] > mincount]
    plottab = plottab_low.append(plottab, ignore_index=True)
    orig_plottab = orig_plottab[orig_plottab['count'] > mincount]
    orig_plottab = orig_plottab_low.append(orig_plottab, ignore_index=True)

    if plot == 'Yes':
        #Plot results
        plt.rcParams.update({'font.size': fontsize})

        if colorlist != 'None':
            colorlist = colorlist
        else:
            colorlist = colorlistX

        #First make figure with smoothed results
        colorcheck = []
        plt.figure(1, figsize=figsize)
        for smp in smplist:
            xdata = plottab['count'][plottab[smp] != 0]
            fdata = plottab[smp][plottab[smp] != 0]

            stype = submeta.loc[smp, var]
            if stype not in colorcheck:
                colorcheck.append(stype)
                plt.plot(xdata, fdata, lw=2, label=stype, color=colorlist[len(colorcheck)-1])
            else:
                plt.plot(xdata, fdata, lw=2, label='_nolegend_', color=colorlist[len(colorcheck) - 1])
        plt.legend(loc='best')
        plt.xlabel('Number of reads')
        plt.xticks(rotation=90)
        plt.ylabel('Diversity')
        if xlim != 'None':
            plt.xlim(xlim[0], xlim[1])
        if ylim != 'None':
            plt.ylim(ylim[0], ylim[1])
        plt.tight_layout()
        if savename != 'None':
            plt.savefig(savename+'_smoothedfig')

        # Then make figure to check if smoothing is ok
        colorcheck = []
        plt.figure(2, figsize=figsize)
        for smp in smplist:
            xdata = plottab['count'][plottab[smp] != 0]
            ydata = orig_plottab[smp][plottab[smp] != 0]
            fdata = plottab[smp][plottab[smp] != 0]

            print(smp, ' ....')
            print('R2= ', r_squared(ydata, fdata))
            print(parameter_dict_low[smp], ' Asymp low= ', parameter_dict_low[smp][0])
            print(parameter_dict[smp], ' Asymp high= ', parameter_dict[smp][1] / parameter_dict[smp][2])
            print('......')

            stype = submeta.loc[smp, var]
            if stype not in colorcheck:
                colorcheck.append(stype)
                plt.plot(xdata, fdata, lw=2, label=stype, color=colorlist[len(colorcheck) - 1])
                plt.plot(xdata, ydata, lw=0, marker='.', label=stype, color=colorlist[len(colorcheck) - 1])
            else:
                plt.plot(xdata, fdata, lw=2, label='_nolegend_', color=colorlist[len(colorcheck) - 1])
                plt.plot(xdata, ydata, lw=0, marker='.', label='_nolegend_', color=colorlist[len(colorcheck) - 1])
        plt.legend(loc='best')
        plt.xlabel('Number of reads')
        plt.xticks(rotation=90)
        plt.ylabel('Diversity')
        if xlim != 'None':
            plt.xlim(xlim[0], xlim[1])
        if ylim != 'None':
            plt.ylim(ylim[0], ylim[1])
        plt.tight_layout()
        if savename != 'None':
            plt.savefig(savename+'_controlfig')

        plt.show()

    else:
        return plottab


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
