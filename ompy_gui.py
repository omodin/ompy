import tkinter as tk
from tkinter.filedialog import askdirectory
from tkinter.filedialog import askopenfilename
from ompy_functions import *
import pandas as pd

def Heatmap(obj, path):
    def run():
        levels = []
        for i in range(len(v)):
            if v[i].get() == 1:
                levels.append(options[i])

        if useLabels.get() == 'Yes':
            labs = True
        else:
            labs = False

        if len(usecbar.get()) > 0 and usecbar.get() != 'None':
            cbar = usecbar.get()
            cbar = cbar.replace(' ', '')
            cbar = cbar.split(',')
            for i in range(len(cbar)):
                cbar[i] = float(cbar[i])
        else:
            cbar = []

        if len(sepCol.get()) > 0 and sepCol.get() != 'None':
            sepc = sepCol.get()
            sepc = sepc.replace(' ', '')
            sepc = sepc.split(',')
            for i in range(len(sepc)):
                sepc[i] = int(sepc[i])
        else:
            sepc = []

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

        plotHeatmap(obj, xAxis=var.get(), levels=levels, subsetLevels=stringlevels, subsetPatterns=stringpatterns,
                        order=order.get(), numberToPlot=number.get(), method='max_sample', nameType=nametype.get(),
                        figSize=(figSizeW.get(), figSizeH.get()), fontSize=fontSize.get(), sepCol=sepc,
                        labels=labs, labelSize=labelSize.get(), cThreshold=ctresh.get(),
                        cMap=colmap.get(), cLinear=colgamma.get(), cBar=cbar, savename=path+'Heatmap')

    def quit():
        master.destroy()

    # Create GUI window
    master = tk.Toplevel()
    master.title('Heatmap')
    master.geometry('600x700')

    # Create scrollbar, root is the frame the all widgets are later placed on
    def onFrameConfigure(canvas):
        ###Reset the scroll region to encompass the inner frame
        canvas.configure(scrollregion=canvas.bbox("all"))
    canvas = tk.Canvas(master, borderwidth=0, background="#bebebe")
    root = tk.Frame(canvas, background="#bebebe")
    vsb = tk.Scrollbar(master, orient="vertical", command=canvas.yview)
    hsb = tk.Scrollbar(master, orient="horizontal", command=canvas.xview)

    canvas.configure(yscrollcommand=vsb.set)
    canvas.configure(xscrollcommand=hsb.set)

    vsb.pack(side="right", fill="y")
    hsb.pack(side="bottom", fill="x")

    canvas.pack(side="left", fill="both", expand=True)
    canvas.create_window((8, 20), window=root, anchor="nw")

    root.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

    ####################
    tk.Label(root, text='Various input options for heatmap').grid(row=0, columnspan=3, sticky=tk.W)
    tk.Label(root, text='-'*120).grid(row=1, columnspan=3, sticky=tk.W)

    # xAxis
    var = tk.StringVar(root, 'None')
    tk.Label(root, text='Enter metadata column for x-axis labels').grid(row=2, columnspan=3, sticky=tk.W)
    tk.Entry(root, textvariable=var).grid(row=3, sticky=tk.W)

    # Input taxonomic levels
    options = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    v = []
    for i in range(len(options)):
        v.append(tk.IntVar())
    tk.Label(root, text='Choose one or two taxonomic levels to include on the y-axis').grid(row=10, columnspan=3, sticky=tk.W)
    tk.Label(root, text='Sequences are grouped based on the lowest taxonomic level').grid(row=11, columnspan=3, sticky=tk.W)
    for val, opt in enumerate(options):
        tk.Checkbutton(root, text=opt, variable=v[val]).grid(row=12+val, sticky=tk.W)

    #Order
    order = tk.StringVar(root, 'None')
    tk.Label(root, text='Specify metadata column used to order the samples on the x-axis').grid(row=20, columnspan=3, sticky=tk.W)
    tk.Entry(root, textvariable=order).grid(row=21, sticky=tk.W)

    #Number to plot
    number = tk.IntVar(root, 20)
    tk.Label(root, text='Specify number of taxa to include in heatmap').grid(row=23, columnspan=3, sticky=tk.W)
    tk.Entry(root, textvariable=number, width=10).grid(row=24, sticky=tk.W)

    #nameType
    nametype = tk.StringVar()
    nametype.set('SV')
    tk.Label(root, text='Specify how unclassified taxa should be named (e.g. SV or OTU)?').grid(row=26, columnspan=3, sticky=tk.W)
    tk.Entry(root, textvariable=nametype, width=10).grid(row=27, sticky=tk.W)

    #Figure dimensions
    tk.Label(root, text='-'*120).grid(row=30, columnspan=3, sticky=tk.W)

    #figSize
    tk.Label(root, text='Specify figure dimensions and text size').grid(row=31, columnspan=3, sticky=tk.W)
    figSizeW = tk.IntVar(root, 14)
    figSizeH = tk.IntVar(root, 10)
    tk.Label(root, text='Width').grid(row=32, sticky=tk.E)
    tk.Entry(root, textvariable=figSizeW).grid(row=32, column=1, sticky=tk.W)
    tk.Label(root, text='Height').grid(row=33, sticky=tk.E)
    tk.Entry(root, textvariable=figSizeH).grid(row=33, column=1, sticky=tk.W)

    #FontSize
    fontSize = tk.IntVar(root, 15)
    tk.Label(root, text='Axis text font size').grid(row=35, sticky=tk.E)
    tk.Entry(root, textvariable=fontSize).grid(row=35, column=1, sticky=tk.W)

    #sepCol
    tk.Label(root, text='-'*120).grid(row=36, columnspan=3, sticky=tk.W)
    sepCol = tk.StringVar()
    tk.Label(root, text='Group samples. Insert numbers of samples after which a blank column should be inserted').grid(row=37, columnspan=3, sticky=tk.W)
    tk.Entry(root, textvariable=sepCol).grid(row=38, sticky=tk.W)
    tk.Label(root, text='(separate values by commas)').grid(row=38, column=1, sticky=tk.W)

    #Data labels
    tk.Label(root, text='-'*120).grid(row=40, columnspan=3, sticky=tk.W)
    tk.Label(root, text='Information about data labels').grid(row=41, sticky=tk.W)

    tk.Label(root, text='Do you want to include data labels in heatmap').grid(row=42, columnspan=3, sticky=tk.W)
    labeloptions = ['Yes', 'No']
    useLabels = tk.StringVar(root, 'Yes')
    for val, opt in enumerate(labeloptions):
        tk.Radiobutton(root, text=opt, variable=useLabels, value=opt).grid(row=43, column=val, sticky=tk.W)

    labelSize = tk.IntVar(root, 12)
    tk.Label(root, text='Label font size').grid(row=45, sticky=tk.W)
    tk.Entry(root, textvariable=labelSize, width=10).grid(row=46, sticky=tk.W)

    ctresh = tk.IntVar(root, 10)
    tk.Label(root, text='Percent relative abundance at which the label text shifts from black to white').grid(row=48, columnspan=3, sticky=tk.W)
    tk.Entry(root, textvariable=ctresh, width=10).grid(row=49, sticky=tk.W)

    #Coloring
    tk.Label(root, text='-'*120).grid(row=50, columnspan=3, sticky=tk.W)
    tk.Label(root, text='Color of heatmap').grid(row=51, columnspan=3, sticky=tk.W)
    colmap = tk.StringVar(root, 'Reds')
    tk.Label(root, text='Colormap').grid(row=52, sticky=tk.E)
    tk.Entry(root, textvariable=colmap).grid(row=52, column=1, sticky=tk.W)
    tk.Label(root, text='(see available colormaps in python)').grid(row=52, column=2, sticky=tk.W)

    colgamma = tk.DoubleVar(root, 0.5)
    tk.Label(root, text='Linearity of colormap').grid(row=55, sticky=tk.E)
    tk.Entry(root, textvariable=colgamma).grid(row=55, column=1, sticky=tk.W)
    tk.Label(root, text='(1=linear change in color)').grid(row=55, column=2, sticky=tk.W)

    tk.Label(root, text='If you want a colorbar showing the scale, specify tick marks on the bar').grid(row=58, columnspan=3, sticky=tk.W)
    usecbar = tk.StringVar(root, 'None')
    tk.Entry(root, textvariable=usecbar).grid(row=59, sticky=tk.W)
    tk.Label(root, text='(the values should be separated by comma)').grid(row=59, column=1, columnspan=2, sticky=tk.W)

    # Input subset based on string patterns
    stringlev = []
    for i in range(len(options)):
        stringlev.append(tk.IntVar())
    tk.Label(root, text='-'*120).grid(row=60, columnspan=3, sticky=tk.W)
    tk.Label(root, text='Subset data based on text patterns').grid(row=61, columnspan=3, sticky=tk.W)
    tk.Label(root, text='Choose which taxonomic levels to search for text').grid(row=62, columnspan=3, sticky=tk.W)
    for val, opt in enumerate(options):
        tk.Checkbutton(root, text=opt, variable=stringlev[val]).grid(row=63+val, sticky=tk.W)

    stringpattern = tk.StringVar()
    tk.Label(root, text='Enter words to search for, separate by comma').grid(row=70, columnspan=3, sticky=tk.W)
    tk.Entry(root, textvariable=stringpattern, width=50).grid(row=71, columnspan=2, sticky=tk.W)
    tk.Label(root, text='-'*120).grid(row=72, columnspan=3, sticky=tk.W)

    # Buttons to run functions
    tk.Button(root, text='Plot heatmap', command=run).grid(row=80)
    tk.Button(root, text='Quit', command=quit).grid(row=80, column=1)
    tk.Label(root, text='-'*120).grid(row=82, columnspan=3, sticky=tk.W)

    root.mainloop()

def Alpha_div(obj, path):
    # Create GUI window
    master = tk.Toplevel()
    master.title('Alpha diversity')
    master.geometry('500x700')

    # Create scrollbar, root is the frame the all widgets are later placed on
    def onFrameConfigure(canvas):
        ###Reset the scroll region to encompass the inner frame
        canvas.configure(scrollregion=canvas.bbox("all"))
    canvas = tk.Canvas(master, borderwidth=0, background="#bebebe")
    root = tk.Frame(canvas, background="#bebebe")
    vsb = tk.Scrollbar(master, orient="vertical", command=canvas.yview)
    hsb = tk.Scrollbar(master, orient="horizontal", command=canvas.xview)

    canvas.configure(yscrollcommand=vsb.set)
    canvas.configure(xscrollcommand=hsb.set)

    vsb.pack(side="right", fill="y")
    hsb.pack(side="bottom", fill="x")

    canvas.pack(side="left", fill="both", expand=True)
    canvas.create_window((8, 20), window=root, anchor="nw")

    root.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))
    ########

    tk.Label(root, text='Show alpha diversity in plots or save as data files').pack(anchor=tk.W)
    tk.Label(root, text='-'*100).pack(anchor=tk.W)

    # Input distmat
    tk.Label(root, text='Are you working with phylogenetic diversity?').pack(anchor=tk.W)
    tk.Label(root, text='If so, select distance matrix file').pack(anchor=tk.W)
    distmat = tk.StringVar(root, 'Select')
    def openDistmat():
        distmat.set(askopenfilename())
    tk.Button(root, textvariable=distmat, command=openDistmat).pack(anchor=tk.W)
    def resetNone():
        distmat.set('Select')
    tk.Button(root, text='Reset selection', command=resetNone).pack(anchor=tk.W)

    # Rarefy
    rare = tk.StringVar()
    rare.set('min')
    tk.Label(root, text='Do you want to rarefy the data?').pack(anchor=tk.W)
    tk.Label(root, text='The default is to rarefy to the smallest sample (min)').pack(anchor=tk.W)
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

        if distmat.get() != 'Select':
            fulldistmat = pd.read_csv(distmat.get(), index_col=0)
            name2save = 'Phyl_alpha_div_fig'
        else:
            fulldistmat = 'None'
            name2save = 'Naive_alpha_div_fig'

        plotDivAlpha(obj, distmat=fulldistmat, rarefy=rarefy, var=var_col.get(), slist='All', order=order.get(), ylog=ylog,
                     colorlist='None', savename=path + name2save)

    tk.Label(root, text='-'*100).pack(anchor=tk.W)
    tk.Label(root, text='The following input is used for plotting figure...').pack(anchor=tk.W)
    tk.Label(root, text='. '*50).pack(anchor=tk.W)

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
            if distmat.get() == 'Select':
                alfa = naiveDivAlpha(obj['tab'], q=q, rarefy=rarefy)
                output[q] = alfa
            else:
                dist = pd.read_csv(distmat.get(), index_col=0)
                alfa = phylDivAlpha(obj['tab'], distmat=dist, q=q, rarefy=rarefy)
                output[q] = alfa

        if distmat.get() == 'Select':
            name2save = 'Naive_alpha_div.csv'
        else:
            name2save = 'Phyl_alpha_div.csv'
        output.to_csv(path + name2save)

    tk.Label(root, text='-'*100).pack(anchor=tk.W)
    tk.Label(root, text='The following input is used to save a csv file with data').pack(anchor=tk.W)
    tk.Label(root, text='. '*50).pack(anchor=tk.W)

    tk.Label(root, text='Specify diversity orders to calculate, use comma to separate numbers').pack(anchor=tk.W)
    qvalues = tk.StringVar()
    tk.Entry(root, textvariable=qvalues).pack(anchor=tk.W)

    # Buttons to save
    tk.Button(root, text='Save alpha diversity data as file', command=run_print).pack(anchor=tk.W)
    tk.Label(root, text='-'*100).pack(anchor=tk.W)

    def quit():
        master.destroy()
    tk.Button(root, text='Quit', command=quit).pack(anchor=tk.W)

def Beta_div(obj, path):

    # Create GUI window
    master = tk.Toplevel()
    master.title('Beta diversity')
    master.geometry('500x700')

    # Create scrollbar, root is the frame the all widgets are later placed on
    def onFrameConfigure(canvas):
        ###Reset the scroll region to encompass the inner frame
        canvas.configure(scrollregion=canvas.bbox("all"))
    canvas = tk.Canvas(master, borderwidth=0, background="#bebebe")
    root = tk.Frame(canvas, background="#bebebe")
    vsb = tk.Scrollbar(master, orient="vertical", command=canvas.yview)
    hsb = tk.Scrollbar(master, orient="horizontal", command=canvas.xview)

    canvas.configure(yscrollcommand=vsb.set)
    canvas.configure(xscrollcommand=hsb.set)

    vsb.pack(side="right", fill="y")
    hsb.pack(side="bottom", fill="x")

    canvas.pack(side="left", fill="both", expand=True)
    canvas.create_window((8, 20), window=root, anchor="nw")

    root.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))
    ######

    tk.Label(root, text='Show beta diversity in PCoA and save pairwise dissimilarities as data files').pack(anchor=tk.W)
    tk.Label(root, text='-------------------------------------------------').pack(anchor=tk.W)

    #Plot or save
    tk.Label(root, text='Do you want to plot PCoA?').pack(anchor=tk.W)
    options = ['Yes', 'No']
    choice = tk.StringVar(root, 'Yes')
    for val, opt in enumerate(options):
        tk.Radiobutton(root, text=opt, variable=choice, value=opt).pack(anchor=tk.W)
    tk.Label(root, text='-------------------------------------------------').pack(anchor=tk.W)

    # Input distmat
    tk.Label(root, text='Are you working with phylogenetic diversity?').pack(anchor=tk.W)
    tk.Label(root, text='If so, select distance matrix file').pack(anchor=tk.W)
    distmat = tk.StringVar(root, 'Select')
    def openDistmat():
        distmat.set(askopenfilename())
    tk.Button(root, textvariable=distmat, command=openDistmat).pack(anchor=tk.W)
    def resetNone():
        distmat.set('Select')
    tk.Button(root, text='Reset selection', command=resetNone).pack(anchor=tk.W)

    # Rarefy
    rare = tk.StringVar()
    rare.set('min')
    tk.Label(root, text='Do you want to rarefy the data?').pack(anchor=tk.W)
    tk.Label(root, text='The default is to rarefy to the smallest sample (min)').pack(anchor=tk.W)
    tk.Label(root, text='You can also specify a number or write None if you do not want to rarefy').pack(anchor=tk.W)
    tk.Entry(root, textvariable=rare).pack(anchor=tk.W)

    # q
    qval = tk.StringVar()
    qval.set(1)
    tk.Label(root, text='Specify diversity order (q)').pack(anchor=tk.W)
    tk.Entry(root, textvariable=qval).pack(anchor=tk.W)

    tk.Label(root, text='-------------------------------------------------').pack(anchor=tk.W)
    tk.Label(root, text='The following input is only used for plotting PCoA').pack(anchor=tk.W)

    var_col = tk.StringVar()
    var_col.set('')
    tk.Label(root, text='Set metadata column heading for color coding of points (required)').pack(anchor=tk.W)
    tk.Entry(root, textvariable=var_col).pack(anchor=tk.W)

    var_marker = tk.StringVar()
    var_marker.set('None')
    tk.Label(root, text='Set metadata column heading for marker type of points').pack(anchor=tk.W)
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
        qsavestring = '_q' + str(q).replace('.', '_')

        if distmat.get() != 'Select':
            dist = pd.read_csv(distmat.get(), index_col=0)
            dis = phylDivBeta(obj['tab'], distmat=dist, q=q, rarefy=rarefy, dis=True)
            name2save = 'Phyl_beta_div'
        else:
            dis = naiveDivBeta(obj['tab'], q=q, rarefy=rarefy, dis=True)
            name2save = 'Naive_beta_div'

        dis.to_csv(path + name2save + qsavestring + '_dissimilarities.csv')
        if choice.get() == 'Yes':
            plotPCoA(dis, obj['meta'], var1=var_col.get(), var2=var_marker.get(), savename=path+name2save+qsavestring+'_PCoA')

    tk.Button(root, text='Calculate beta', command=dis_func).pack(anchor=tk.W)

    def quit():
        master.destroy()
    tk.Button(root, text='Quit', command=quit).pack(anchor=tk.W)

def Calc_phyl_dist(obj, path):
    # Create GUI window
    master = tk.Toplevel()
    master.title('Phylogenetic distances')
    master.geometry('500x400')

    # Create scrollbar, root is the frame the all widgets are later placed on
    def onFrameConfigure(canvas):
        ###Reset the scroll region to encompass the inner frame
        canvas.configure(scrollregion=canvas.bbox("all"))
    canvas = tk.Canvas(master, borderwidth=0, background="#bebebe")
    root = tk.Frame(canvas, background="#bebebe")
    vsb = tk.Scrollbar(master, orient="vertical", command=canvas.yview)
    hsb = tk.Scrollbar(master, orient="horizontal", command=canvas.xview)

    canvas.configure(yscrollcommand=vsb.set)
    canvas.configure(xscrollcommand=hsb.set)

    vsb.pack(side="right", fill="y")
    hsb.pack(side="bottom", fill="x")

    canvas.pack(side="left", fill="both", expand=True)
    canvas.create_window((8, 20), window=root, anchor="nw")

    root.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))
    ###########

    # Name of distmat
    distmat_name = tk.StringVar()
    distmat_name.set('phyl_dist')
    tk.Label(root, text='Specify name of the file to be saved (e.g. phyl_dist').pack(anchor=tk.W)
    tk.Label(root, text='(No need to add .csv, it will be added automatically)').pack(anchor=tk.W)
    tk.Entry(root, textvariable=distmat_name).pack(anchor=tk.W)

    def save_func():
        savename = distmat_name.get()
        phylDistMat(obj['seq'], savename=path+savename)

    tk.Label(root, text='The calculation may take quite long time').pack(anchor=tk.W)
    tk.Button(root, text='Calculate and save file', command=save_func).pack(anchor=tk.W)
    tk.Label(root, text='-'*70).pack(anchor=tk.W)

    def quit():
        master.destroy()
    tk.Button(root, text='Quit', command=quit).pack(anchor=tk.W)

def Subsetting(obj, path):

    # Create GUI window
    master = tk.Toplevel()
    master.title('Beta diversity')
    master.geometry('500x600')

    # Create scrollbar, root is the frame the all widgets are later placed on
    def onFrameConfigure(canvas):
        ###Reset the scroll region to encompass the inner frame
        canvas.configure(scrollregion=canvas.bbox("all"))
    canvas = tk.Canvas(master, borderwidth=0, background="#bebebe")
    root = tk.Frame(canvas, background="#bebebe")
    vsb = tk.Scrollbar(master, orient="vertical", command=canvas.yview)
    hsb = tk.Scrollbar(master, orient="horizontal", command=canvas.xview)

    canvas.configure(yscrollcommand=vsb.set)
    canvas.configure(xscrollcommand=hsb.set)

    vsb.pack(side="right", fill="y")
    hsb.pack(side="bottom", fill="x")

    canvas.pack(side="left", fill="both", expand=True)
    canvas.create_window((8, 20), window=root, anchor="nw")

    root.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))
    ######

    # Top SVs
    top_svs = tk.IntVar()
    tk.Label(root, text='Subset to most abundant sequences').pack(anchor=tk.W)
    tk.Label(root, text='Write number of sequences to keep').pack(anchor=tk.W)
    tk.Entry(root, textvariable=top_svs).pack(anchor=tk.W)

    def subset_top_svs():
        obj_sub = subsetTopSVs(obj, top_svs.get())
        returnFiles(obj_sub, path=path, sep='\t')
    tk.Button(root, text='Subset and save files', command=subset_top_svs).pack(anchor=tk.W)

    tk.Label(root, text='-'*100).pack(anchor=tk.W)
    # Merge samples
    meta_h = tk.StringVar()
    tk.Label(root, text='Merge sample based on metadata column heading').pack(anchor=tk.W)
    tk.Label(root, text='Write column heading').pack(anchor=tk.W)
    tk.Entry(root, textvariable=meta_h).pack(anchor=tk.W)

    def merge_smps():
        obj_sub = mergeSamples(obj, var=meta_h.get())
        returnFiles(obj_sub, path=path, sep='\t')
    tk.Button(root, text='Merge samples and save files', command=merge_smps).pack(anchor=tk.W)

    def quit():
        master.destroy()
    tk.Button(root, text='Quit', command=quit).pack(anchor=tk.W)

def startwindow():

    # Functions that specify what to do with input choices
    def choose():
        if sep_name.get() == 'tab':
            separator = '\t'
        else:
            separator = sep_name.get()

        obj = loadFiles(path='', tab=table_name.get(), fasta=seq_name.get(), meta=meta_name.get(), sep=separator)

        if v.get() == 'Heatmap':
            Heatmap(obj, path_name.get()+'/')
        elif v.get() == 'Alpha_div':
            Alpha_div(obj, path_name.get()+'/')
        elif v.get() == 'Beta_div':
            Beta_div(obj, path_name.get()+'/')
        elif v.get() == 'Calculate_phyl_dist':
            Calc_phyl_dist(obj, path_name.get()+'/')
        elif v.get() == 'Subset_data':
            Subsetting(obj, path_name.get()+'/')

    def quit():
        master.destroy()

    # Start GUI window
    master = tk.Tk()
    master.title('Start window')
    master.geometry('500x600')

    # Create scrollbar, root is the frame the all widgets are later placed on
    def onFrameConfigure(canvas):
        ###Reset the scroll region to encompass the inner frame
        canvas.configure(scrollregion=canvas.bbox("all"))
    canvas = tk.Canvas(master, borderwidth=0, background="#bebebe")
    root = tk.Frame(canvas, background="#bebebe")
    vsb = tk.Scrollbar(master, orient="vertical", command=canvas.yview)
    hsb = tk.Scrollbar(master, orient="horizontal", command=canvas.xview)

    canvas.configure(yscrollcommand=vsb.set)
    canvas.configure(xscrollcommand=hsb.set)

    vsb.pack(side="right", fill="y")
    hsb.pack(side="bottom", fill="x")

    canvas.pack(side="left", fill="both", expand=True)
    canvas.create_window((8, 20), window=root, anchor="nw")

    root.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

    ##########
    tk.Label(root, text='Select input files').grid(row=1, sticky=tk.W)

    table_name = tk.StringVar(root, ' ')
    def openTable():
        table_name.set(askopenfilename())
    tk.Button(root, text='Frequency table', command=openTable).grid(row=2, sticky=tk.W)
    tk.Label(root, textvariable=table_name).grid(row=2, column=1, sticky=tk.W)

    seq_name = tk.StringVar(root, ' ')
    def openSeq():
        seq_name.set(askopenfilename())
    tk.Button(root, text='Fasta file', command=openSeq).grid(row=3, sticky=tk.W)
    tk.Label(root, textvariable=seq_name).grid(row=3, column=1, sticky=tk.W)

    meta_name = tk.StringVar(root, ' ')
    def openMeta():
        meta_name.set(askopenfilename())
    tk.Button(root, text='Meta data', command=openMeta).grid(row=4, sticky=tk.W)
    tk.Label(root, textvariable=meta_name).grid(row=4, column=1, sticky=tk.W)

    tk.Label(root, text='Which type of separator was used in the table and meta files?').grid(row=5, columnspan=2, sticky=tk.W)
    sep_name = tk.StringVar(root, 'tab')
    optionsSep = ['tab', ',', ';']
    for val, opt in enumerate(optionsSep):
        tk.Radiobutton(root, text=opt, variable=sep_name, value=opt).grid(row=10+val, sticky=tk.W)

    #Specify path for output files
    path_name = tk.StringVar(root, ' ')
    def openFolder():
        path_name.set(askdirectory())
    tk.Button(root, text='Folder for output files', command=openFolder).grid(row=16, sticky=tk.W)
    tk.Label(root, textvariable=path_name).grid(row=16, column=1, sticky=tk.W)

    # Choices of analysis
    tk.Label(root, text='-'*100).grid(row=30, columnspan=2, sticky=tk.W)
    tk.Label(root, text='Choose a task').grid(row=31, sticky=tk.W)

    v = tk.StringVar()
    options = ['Heatmap', 'Alpha_div', 'Beta_div', 'Calculate_phyl_dist', 'Subset_data']
    for val, opt in enumerate(options):
        tk.Radiobutton(root, text=opt, variable=v, value=opt).grid(row=40+val, sticky=tk.W)

    # Buttons that connects to functions
    tk.Label(root, text='-'*100).grid(row=51, columnspan=2, sticky=tk.W)
    tk.Button(root, text='Choose', command=choose).grid(row=52, sticky=tk.W)
    tk.Button(root, text='Quit', command=quit).grid(row=52, column=1, sticky=tk.W)


    root.mainloop()

startwindow()


