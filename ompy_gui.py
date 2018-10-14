import tkinter as tk
from ompy_functions import *
import pandas as pd

def Heatmap(obj):
    def run():
        levels = []
        for i in range(len(v)):
            if v[i].get() == 1:
                levels.append(options[i])

        stringlevels = []
        for i in range(len(stringlev)):
            if stringlev[i].get() == 1:
                stringlevels.append(options[i])
        if len(stringlevels) <1:
            stringlevels = 'None'

        if stringpattern.get() != 'None':
            stringpatterns = stringpattern.get().replace(' ', '')
            stringpatterns = stringpatterns.split(',')
        else:
            stringpatterns = 'None'

        plotHeatmap(obj, var=var.get(), levels=levels, subsetLevels=stringlevels, subsetPatterns=stringpatterns, order=order.get(),
                    numberToPlot=int(number.get()), nameType=nametype.get(), labelSize=int(labelSize.get()), fontSize=int(fontSize.get()))

    def quit():
        root.destroy()

    # Create GUI window
    root = tk.Toplevel()
    root.title('Heatmap')

    # Input var
    var = tk.StringVar()
    var.set('None')
    tk.Label(root, text='Enter metadata column for x-axis labels').pack(anchor=tk.W)
    tk.Entry(root, textvariable=var).pack(anchor=tk.W)

    # Input taxonomic levels
    options = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    v = []
    for i in range(len(options)):
        v.append(tk.IntVar())
    tk.Label(root, text='Choose one or two taxonomic levels to include on the y-axis').pack(anchor=tk.W)
    tk.Label(root, text='sequences are grouped based on the lowest').pack(anchor=tk.W)
    for val, opt in enumerate(options):
        tk.Checkbutton(root, text=opt, variable=v[val]).pack(anchor=tk.W)

    #Order
    order = tk.StringVar()
    order.set('None')
    tk.Label(root, text='Specify metadata column used to order the samples on the x-axis').pack(anchor=tk.W)
    tk.Entry(root, textvariable=order).pack(anchor=tk.W)

    #Number to plot
    number = tk.StringVar()
    number.set('20')
    tk.Label(root, text='Specify number of taxa to include in heatmap').pack(anchor=tk.W)
    tk.Entry(root, textvariable=number).pack(anchor=tk.W)

    #nameType
    nametype = tk.StringVar()
    nametype.set('SV')
    tk.Label(root, text='How to name unclassified taxa (e.g. SV or OTU)?').pack(anchor=tk.W)
    tk.Entry(root, textvariable=nametype).pack(anchor=tk.W)

    #LabelSize
    labelSize = tk.StringVar()
    labelSize.set('1')
    tk.Label(root, text='Specify size of data labels in heatmap').pack(anchor=tk.W)
    tk.Entry(root, textvariable=labelSize).pack(anchor=tk.W)

    #FontSize
    fontSize = tk.StringVar()
    fontSize.set('15')
    tk.Label(root, text='Specify size of text on axes').pack(anchor=tk.W)
    tk.Entry(root, textvariable=fontSize).pack(anchor=tk.W)

    # Input subset based on string patterns
    stringlev = []
    for i in range(len(options)):
        stringlev.append(tk.IntVar())
    tk.Label(root, text='-------------------------------').pack(anchor=tk.W)
    tk.Label(root, text='Subset data based on text patters').pack(anchor=tk.W)
    tk.Label(root, text='Search for text within taxonomic levels').pack(anchor=tk.W)
    for val, opt in enumerate(options):
        tk.Checkbutton(root, text=opt, variable=stringlev[val]).pack(anchor=tk.W)

    stringpattern = tk.StringVar()
    tk.Label(root, text='Enter words to search for, separate by comma').pack(anchor=tk.W)
    tk.Entry(root, textvariable=stringpattern).pack(anchor=tk.W)
    tk.Label(root, text='-------------------------------').pack(anchor=tk.W)

    # Buttons to run functions
    tk.Button(root, text='Plot heatmap', command=run).pack(anchor=tk.W)
    tk.Button(root, text='Quit', command=quit).pack(anchor=tk.W)


    root.mainloop()

def Alpha_div(obj, path):
    # Create GUI window
    root = tk.Toplevel()
    root.title('Alpha diversity')
    tk.Label(root, text='Show alpha diversity in plots or save as data files').pack(anchor=tk.W)

    # Input distmat
    distmat = tk.StringVar()
    distmat.set('None')
    tk.Label(root, text='Are you working with phylogenetic diversity?').pack(anchor=tk.W)
    tk.Label(root, text='If so, specify name of distance matrix file').pack(anchor=tk.W)
    tk.Entry(root, textvariable=distmat).pack(anchor=tk.W)

    # Rarefy
    rare = tk.StringVar()
    rare.set('min')
    tk.Label(root, text='Do you want to rarefy the data?').pack(anchor=tk.W)
    tk.Label(root, text='The default is to rarefy to the smallest sample').pack(anchor=tk.W)
    tk.Label(root, text='You can also specify a number or write None if you do not want to rarefy').pack(anchor=tk.W)
    tk.Entry(root, textvariable=rare).pack(anchor=tk.W)

    #Plotting
    def run_plot():
        if rare.get() == 'None' or rare.get() == 'min':
            rarefy = rare.get()
        else:
            rarefy = int(rare.get())

        if y_v.get() == 'Yes':
            ylog = True
        else:
            ylog = False

        plotDivAlpha(obj, distmat=distmat.get(), rarefy=rarefy, var=var_col.get(), slist='All', order=order.get(), ylog=ylog,
                     colorlist='None', savename='None')

    tk.Label(root, text='------------------------').pack(anchor=tk.W)
    tk.Label(root, text='The following input is used for plotting figure...').pack(anchor=tk.W)

    # var to use for color coding
    var_col = tk.StringVar()
    var_col.set('None')
    tk.Label(root, text='Specify metadata column heading to use for color coding').pack(anchor=tk.W)
    tk.Entry(root, textvariable=var_col).pack(anchor=tk.W)

    #Order
    order = tk.StringVar()
    order.set('None')
    tk.Label(root, text='Specify metadata column used to order the samples on the x-axis').pack(anchor=tk.W)
    tk.Entry(root, textvariable=order).pack(anchor=tk.W)

    #Semi log y-axis
    options = ['Yes', 'No']
    tk.Label(root, text='Use logarithmic y-axis?').pack(anchor=tk.W)
    y_v = tk.StringVar()
    for opt in options:
        tk.Radiobutton(root, text=opt, variable=y_v, value=opt).pack(anchor=tk.W)

    # Buttons to plot
    tk.Button(root, text='Plot alpha diversity figures', command=run_plot).pack(anchor=tk.W)

    ## Printing
    def run_print():
        qlist = qvalues.get().replace(' ', '')
        qlist = qlist.split(',')
        qnumbers = []
        for q in qlist:
            qnumbers.append(float(q))

        if rare.get() == 'None' or rare.get() == 'min':
            rarefy = rare.get()
        else:
            rarefy = int(rare.get())

        output = pd.DataFrame(0, index=obj['tab'].columns, columns=qnumbers)
        for q in qnumbers:
            if distmat.get() == 'None':
                alfa = naiveDivAlpha(obj['tab'], q=q, rarefy=rarefy)
                output[q] = alfa
            else:
                alfa = phylDivAlpha(obj['tab'], distmat=distmat, q=q, rarefy=rarefy)
                output[q] = alfa

        output.to_csv(path+'Alpha_diversity.csv')

    tk.Label(root, text='-------------------------').pack(anchor=tk.W)
    tk.Label(root, text='Save data file').pack(anchor=tk.W)
    tk.Label(root, text='Specify diversity orders to print').pack(anchor=tk.W)
    tk.Label(root, text='use comma to separate numbers').pack(anchor=tk.W)
    qvalues = tk.StringVar()
    tk.Entry(root, textvariable=qvalues).pack(anchor=tk.W)

    # Buttons to save
    tk.Button(root, text='Save alpha diversity data as file', command=run_print).pack(anchor=tk.W)

    def quit():
        root.destroy()
    tk.Button(root, text='Quit', command=quit).pack(anchor=tk.W)

def Beta_div(obj, path):

    # Create GUI window
    root = tk.Toplevel()
    root.title('Beta diversity')
    tk.Label(root, text='Show beta diversity in PCoA and save pairwise dissimilarities as data files').pack(anchor=tk.W)
    tk.Label(root, text='-------------------------------------------------').pack(anchor=tk.W)

    #Plot or save
    options = ['Plot_PCoA', 'Save_file']
    tk.Label(root, text='What do you want to do? Plot, save, or both?').pack(anchor=tk.W)
    choice = []
    for val, opt in enumerate(options):
        choice.append(tk.IntVar())
        choice[val].set(1)
        tk.Checkbutton(root, text=opt, variable=choice[val]).pack(anchor=tk.W)
    tk.Label(root, text='-------------------------------------------------').pack(anchor=tk.W)

    # Input distmat
    distmat = tk.StringVar()
    distmat.set('None')
    tk.Label(root, text='Are you working with phylogenetic diversity?').pack(anchor=tk.W)
    tk.Label(root, text='If so, specify name of distance matrix file').pack(anchor=tk.W)
    tk.Entry(root, textvariable=distmat).pack(anchor=tk.W)

    # Rarefy
    rare = tk.StringVar()
    rare.set('min')
    tk.Label(root, text='Do you want to rarefy the data?').pack(anchor=tk.W)
    tk.Label(root, text='The default is to rarefy to the smallest sample').pack(anchor=tk.W)
    tk.Label(root, text='You can also specify a number or write None if you do not want to rarefy').pack(anchor=tk.W)
    tk.Entry(root, textvariable=rare).pack(anchor=tk.W)

    # q
    qval = tk.StringVar()
    qval.set(1)
    tk.Label(root, text='Specify diversity order (q)').pack(anchor=tk.W)
    tk.Entry(root, textvariable=qval).pack(anchor=tk.W)

    tk.Label(root, text='-------------------------------------------------').pack(anchor=tk.W)
    tk.Label(root, text='Input required for plotting PCoA').pack(anchor=tk.W)

    var_col = tk.StringVar()
    var_col.set('')
    tk.Label(root, text='Set metadata column heading for color coding of points (required)').pack(anchor=tk.W)
    tk.Entry(root, textvariable=var_col).pack(anchor=tk.W)

    var_marker = tk.StringVar()
    var_marker.set('None')
    tk.Label(root, text='Set metadata column heading for marker type of points (optional)').pack(anchor=tk.W)
    tk.Entry(root, textvariable=var_marker).pack(anchor=tk.W)

    order = tk.StringVar()
    order.set('None')
    tk.Label(root, text='Specify metadata column used to order the samples in the legend').pack(anchor=tk.W)
    tk.Entry(root, textvariable=order).pack(anchor=tk.W)

    def dis_func():
        if rare.get() == 'None' or rare.get() == 'min':
            rarefy = rare.get()
        else:
            rarefy = int(rare.get())

        q = float(qval.get())

        if distmat.get() == 'None':
            dis = naiveDivBeta(obj['tab'], q=q, rarefy=rarefy, dis=True)
        else:
            dist = pd.read_csv(path+distmat.get(), index_col=0)
            dis = phylDivBeta(obj['tab'], distmat=dist, q=q, rarefy=rarefy, dis=True)
        if choice[0].get() == 1:
            plotPCoA(dis, obj['meta'], var1=var_col.get(), var2=var_marker.get())
        if choice[1].get() == 1:
            dis.to_csv(path+'dissimilarities.csv')

    tk.Button(root, text='Run', command=dis_func).pack(anchor=tk.W)

    def quit():
        root.destroy()
    tk.Button(root, text='Quit', command=quit).pack(anchor=tk.W)

def Calc_phyl_dist(obj, path):
    # Create GUI window
    root = tk.Toplevel()
    root.title('Phylogenetic distances')

    # Name of distmat
    distmat_name = tk.StringVar()
    distmat_name.set('phyl_dist')
    tk.Label(root, text='Specify name of file with pairwise distances between sequences').pack(anchor=tk.W)
    tk.Label(root, text='(no need to add .csv, it was be added automatically)').pack(anchor=tk.W)
    tk.Entry(root, textvariable=distmat_name).pack(anchor=tk.W)

    def save_func():
        savename = distmat_name.get()
        phylDistMat(obj['seq'], 'csv', savename=path+savename)

    tk.Button(root, text='Calculate and save file', command=save_func).pack(anchor=tk.W)

def startwindow():

    # Functions that specify what to do with input choices
    def choose():
        if sep_name.get() == 'tab':
            separator = '\t'
        else:
            separator = sep_name.get()
        obj = loadFiles(path=path_name.get(), tab=table_name.get(), fasta=fasta_name.get(), meta=meta_name.get(), sep=separator)
        if v.get() == 'Heatmap':
            Heatmap(obj)
        elif v.get() == 'Alpha_div':
            Alpha_div(obj, path_name.get())
        elif v.get() == 'Beta_div':
            Beta_div(obj, path_name.get())
        elif v.get() == 'Calculate_phyl_dist':
            Calc_phyl_dist(obj, path_name.get())

    def quit():
        root.destroy()

    # Start GUI window
    root = tk.Tk()
    root.title('Start window')

    # Prepare input for loadFiles function
    path_name = tk.StringVar(root, '/home/om/ompy_input/')
    tk.Label(root, text='Path to input files').pack(anchor=tk.W)
    tk.Entry(root, textvariable=path_name).pack(anchor=tk.W)

    table_name = tk.StringVar()
    table_name.set('GUItable.csv')
    tk.Label(root, text='Input name of frequency table file').pack(anchor=tk.W)
    tk.Entry(root, textvariable=table_name).pack(anchor=tk.W)

    fasta_name = tk.StringVar()
    fasta_name.set('GUIseq.fa')
    tk.Label(root, text='Input name of fasta file').pack(anchor=tk.W)
    tk.Entry(root, textvariable=fasta_name).pack(anchor=tk.W)

    meta_name = tk.StringVar()
    meta_name.set('GUImetadata.csv')
    tk.Label(root, text='Input name of metadata file').pack(anchor=tk.W)
    tk.Entry(root, textvariable=meta_name).pack(anchor=tk.W)

    sep_name = tk.StringVar()
    sep_name.set('tab')
    tk.Label(root, text='Input separator type (, ; tab)').pack(anchor=tk.W)
    tk.Entry(root, textvariable=sep_name).pack(anchor=tk.W)

    # Choices of analysis
    tk.Label(root, text='----------------').pack(anchor=tk.W)
    tk.Label(root, text='Choose a task').pack(anchor=tk.W)

    v = tk.StringVar()
    options = ['Heatmap', 'Alpha_div', 'Beta_div', 'Calculate_phyl_dist']

    for opt in options:
        tk.Radiobutton(root, text=opt, variable=v, value=opt).pack(anchor=tk.W)

    # Buttons that connects to functions
    tk.Button(root, text='Choose', command=choose).pack(anchor=tk.W)
    tk.Button(root, text='Quit', command=quit).pack(anchor=tk.W)

    root.mainloop()

startwindow()

