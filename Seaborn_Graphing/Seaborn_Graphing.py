import seaborn as sns
import matplotlib as mpl
import numpy as np
from decimal import *
import DCG_Utilities as dcgutils
import scipy
## Allows for easy text manipulation in illustrator by changing this
mpl.rcParams['pdf.fonttype'] = 42

import matplotlib.pyplot as plt
import pandas as pd

font = {'family' : 'Arial',
        'weight' : 'bold',
        'size' : 14}
mpl.rc('font',**font)


## Setup figure Style
sns.set(context="paper", font_scale=2, style = "white")
sns.despine()

class PlottingFunctions:
    def __init__(self, LFQData, CurveData, LabelType, **kwargs):
        self._LFQDataFrame = pd.read_table(LFQData)
        self._CurveDataFrame = pd.read_table(CurveData)
        self._LabelType = LabelType
        return super().__init__(**kwargs)

    def GenerateAllFigures(self):
        GenerateFigure("PAG1_Minus_Col0_Minus.txt", "PAG1_Minus_Col0_Minus_Curve.txt", "PAG1_Minus_Col0_Minus","CPPips","kde", True)
        GenerateFigure("PAG1_Plus_Col0_Plus.txt", "PAG1_Plus_Col0_Plus_Curve.txt", "PAG1_Plus_Col0_Plus","CPPips","kde", True)
        GenerateFigure("RPT4a_Minus_Col0_Minus.txt", "RPT4a_Minus_Col0_Minus_Curve.txt", "RPT4a_Minus_Col0_Minus","RPPips","kde", True)
        GenerateFigure("RPT4a_Plus_Col0_Plus.txt", "RPT4a_Plus_Col0_Plus_Curve.txt", "RPT4a_Plus_Col0_Plus","RPPips","kde", True)
        #GenerateFigure("RPT4b_Minus_Col0_Minus.txt", "RPT4b_Minus_Col0_Minus_curve.txt", "RPT4b_Minus_Col0_Minus","RPPips","kde", True)
        #GenerateFigure("RPT4b_Plus_Col0_Plus.txt", "RPT4b_Plus_Col0_Plus_curve.txt", "RPT4b_Plus_Col0_Plus","RPPips","kde", True)
        test = self
    
    def GenerateFigure(self, outputFileName, plotstyle="kde", shouldLabel = False):
        ##Read in Data
        df = self._LFQDataFrame
        ##Read in Curve
        df2 = self._CurveDataFrame
        ## NOTES
        ## replace NaN in subcomplex labels with unlabeled so that it plots correctly, and then setup color
        ## hue so that it turns grey
        ## 
        df = df.replace({'SubComplexAnnotation': { np.NaN : 'Unlabeled'}})


        ## Think about doing a 1D density plot on the axes instead of a contour plot

        # Subset the Volcano Plot by ProteasomeSubComplex

        RP = df.loc[df['SubComplexAnnotation'] == 'RP']
        CP = df.loc[df['SubComplexAnnotation'] == 'CP']
        PIPS_CP = df.loc[df['SubComplexAnnotation'] == 'PIP_CP']
        PIPS_RP = df.loc[df['SubComplexAnnotation'] == 'PIP_RP']
        PIP = df.loc[df['SubComplexAnnotation'] == 'PIP']

        ## Setup XKCD Colors
        colors = ["cool grey", "windows blue", "pale red", "medium green", "orange", "black"]
        ## Parse Colors
        colors = sns.xkcd_palette(colors);

        ##SetupCusotmPalleteDictionary
        ##if self._LabelType=="AllPips":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[3], PIP_RP=colors[3])
        ##if self._LabelType=="CPPips":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[3], PIP_RP=colors[0])
        ##if self._LabelType=="RPPips":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[0], PIP_RP=colors[3])
        ##if self._LabelType=="ECM29":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[0], PIP_RP=colors[0])
        ##if self._LabelType=="CPRP":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[0], PIP_CP=colors[0], PIP_RP=colors[0])
    
        ##Rick's Palette
        custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[5], PIP_CP=colors[3], PIP_RP=colors[4])


        ## Setup Seaborn facetgrid with lmplot (no line)
        facetgrid = sns.lmplot('Difference','-Log(P-value)', data=df, hue='SubComplexAnnotation', fit_reg=False, palette = custompalette, scatter_kws={'s' : 100}, legend=False, size = 10, aspect = 1)

        ## Force draw of points from facetgrid by axcessing facetgrids ax, and then setting collection zorder
        plt.setp(facetgrid.ax.collections, zorder =100)

        ## Setup Nice Color Pallete for density plots
        red = sns.color_palette("Reds")[-2]
        blue = sns.color_palette("Blues")[-2]
       
        ## Add 2D Density Plot for CP
        sns.kdeplot(CP['Difference'], CP['-Log(P-value)'], cmap="Reds", bw=0.5, shade=True, shade_lowest=False, alpha = 0.5, ax=facetgrid.ax)
        
        ##Add Difference Density Plot for CP if kde style specified
        if plotstyle=="kde":
            sns.kdeplot(CP['Difference'], label='CP Density', color =red, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5)
        
        ##Add 2D Density Plot for RP
        sns.kdeplot(RP['Difference'], RP['-Log(P-value)'], cmap="Blues", bw=0.5, shade=True, shade_lowest=False, alpha = 0.5, ax=facetgrid.ax)
        
        ##Add Difference Density Plot for RP if kde style specified
        if plotstyle=="kde":
            sns.kdeplot(RP['Difference'], label='RP Density', color =blue, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5)

        ## Use curve data to plot significance curve
        plt.plot(df2.iloc[:,0],df2.iloc[:,1], color = "black", label = "SAM Significance")
    
        ## Setup Legend
        facetgrid.fig.get_axes()[0].legend(loc='upper left', title = 'Key')

        ## Setup Axes and axes limit
        myaxes = facetgrid.axes
        myaxes[0,0].set_xlim(-8,13)
        myaxes[0,0].set_ylim(0,7.5) 

        ## Label Pips if needed
        if shouldLabel:
            if labelstyle=="CPPips":
                self.LabelFigure(PIPS_CP, facetgrid.ax, 0.1, colors)
                self.LabelFigure(PIP, facetgrid.ax, 0.1, colors)
            if labelstyle=="RPPips":
                self.LabelFigure(PIPS_RP, facetgrid.ax, 0.1, colors)
                self.LabelFigure(PIP, facetgrid.ax, 0.1, colors)
            if labelstyle=="AllPips":
                self.LabelFigure(PIPS_CP, facetgrid.ax, 0.1, colors)
                self.LabelFigure(PIPS_RP, facetgrid.ax, 0.1, colors)
                self.LabelFigure(PIP, facetgrid.ax, 0.1, colors)
            if labelstyle=="ECM29":
                self.LabelFigure(PIP, facetgrid.ax, 0.1, colors)
        ## Save Figure
        facetgrid.savefig(outputFileName + ".pdf");
        facetgrid.savefig(outputFileName + ".png");

    def GenerateRPT4abFigure(self, outputFileName, plotstyle="kde", shouldLabel = False):
        ##Read in Data
        df = self._LFQDataFrame
        ##Read in Curve
        df2 = self._CurveDataFrame
        ## NOTES
        ## replace NaN in subcomplex labels with unlabeled so that it plots correctly, and then setup color
        ## hue so that it turns grey
        ##

        df = df.replace({'SubComplexAnnotation': { np.NaN : 'Unlabeled'}})

        ##SubComplexAnnotation
        ##If its not RP,CP,PIP_CP, PIP_RP, or PIP, then label Unlabeled
        ##If its RPT4a, label RPT4a
        ##If its RPT4b, label RPT4b


        ## Think about doing a 1D density plot on the axes instead of a contour plot

        # Subset the Volcano Plot by ProteasomeSubComplex

        RP = df.loc[df['SubComplexAnnotation'] == 'RP']
        CP = df.loc[df['SubComplexAnnotation'] == 'CP']
        PIPS_CP = df.loc[df['SubComplexAnnotation'] == 'PIP_CP']
        PIPS_RP = df.loc[df['SubComplexAnnotation'] == 'PIP_RP']
        PIP = df.loc[df['SubComplexAnnotation'] == 'PIP']

        ## Setup XKCD Colors
        colors = ["cool grey", "windows blue", "pale red", "medium green"]
        ## Parse Colors
        colors = sns.xkcd_palette(colors);

        ##SetupCusotmPalleteDictionary
        if self._LabelType=="AllPips":
            custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[3], PIP_RP=colors[3])
        if self._LabelType=="CPPips":
            custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[3], PIP_RP=colors[0])
        if self._LabelType=="RPPips":
            custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[0], PIP_RP=colors[3])
        if self._LabelType=="ECM29":
            custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[0], PIP_RP=colors[0])
    
    



        ## Setup Seaborn facetgrid with lmplot (no line)
        facetgrid = sns.lmplot('Difference','-Log(P-value)', data=df, hue='SubComplexAnnotation', fit_reg=False, palette = custompalette, scatter_kws={'s' : 100}, legend=False, size = 10, aspect = 2)

        ## Force draw of points from facetgrid by axcessing facetgrids ax, and then setting collection zorder
        plt.setp(facetgrid.ax.collections, zorder =100)

        ## Setup Nice Color Pallete for density plots
        red = sns.color_palette("Reds")[-2]
        blue = sns.color_palette("Blues")[-2]
       
        ## Add 2D Density Plot for CP
        sns.kdeplot(CP['Difference'], CP['-Log(P-value)'], cmap="Reds", bw=0.5, shade=True, shade_lowest=False, alpha = 0.5, ax=facetgrid.ax)
        
        ##Add Difference Density Plot for CP if kde style specified
        if plotstyle=="kde":
            sns.kdeplot(CP['Difference'], label='CP Density', color =red, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5)
        
        ##Add 2D Density Plot for RP
        sns.kdeplot(RP['Difference'], RP['-Log(P-value)'], cmap="Blues", bw=0.5, shade=True, shade_lowest=False, alpha = 0.5, ax=facetgrid.ax)
        
        ##Add Difference Density Plot for RP if kde style specified
        if plotstyle=="kde":
            sns.kdeplot(RP['Difference'], label='RP Density', color =blue, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5)

        ## Use curve data to plot significance curve
        plt.plot(df2.iloc[:,0],df2.iloc[:,1], color = "black", label = "SAM Significance")
    
        ## Setup Legend
        facetgrid.fig.get_axes()[0].legend(loc='upper left', title = 'Key')

        ## Setup Axes and axes limit
        myaxes = facetgrid.axes
        myaxes[0,0].set_xlim(-8,13)
        myaxes[0,0].set_ylim(0,7.5)

        ## Label Pips if needed
        if shouldLabel:
            if labelstyle=="CPPips":
                self.LabelFigure(PIPS_CP, facetgrid.ax, 0.1, colors)
                self.LabelFigure(PIP, facetgrid.ax, 0.1, colors)
            if labelstyle=="RPPips":
                self.LabelFigure(PIPS_RP, facetgrid.ax, 0.1, colors)
                self.LabelFigure(PIP, facetgrid.ax, 0.1, colors)
            if labelstyle=="AllPips":
                self.LabelFigure(PIPS_CP, facetgrid.ax, 0.1, colors)
                self.LabelFigure(PIPS_RP, facetgrid.ax, 0.1, colors)
                self.LabelFigure(PIP, facetgrid.ax, 0.1, colors)
            if labelstyle=="ECM29":
                self.LabelFigure(PIP, facetgrid.ax, 0.1, colors)
        ## Save Figure
        facetgrid.savefig(outputFileName + ".pdf");
        facetgrid.savefig(outputFileName + ".png");

    def GetSignificantInteractorsThatAreNotProteasomeSubunits(self, outputfilename):
        dftoProc = self._LFQDataFrame 
        ##Logic for not an interactor
        ## Significant as in +
        ## >0 in difference (works as long as volcano plot is setup correctly
        ## Doesn't Have a ProteasomeID
        #print(dftoProc)
        #Interactors = dftoProc.loc[dftoProc['ProteasomeAnnotation'].isnull()]
        Interactors = dftoProc.loc[dftoProc['SubComplexAnnotation'] != 'CP']
        Interactors = Interactors.loc[Interactors['SubComplexAnnotation'] != 'RP']
        if ('+' in Interactors['Significant'].values):
            Interactors = Interactors.loc[Interactors['Significant'] == '+']
            Interactors = Interactors.loc[Interactors['Difference'] >= 0]
            Outputlist = ['Fasta headers','Difference', '-Log(P-value)', 'ProteasomeAnnotation']
            #print(Interactors)
            dataToWrite = dcgutils.Utilities.FilterDataFramebyInclusionList(Interactors, Outputlist, sortbycolumnname = 'Difference')
            dataToWrite = dataToWrite.reset_index()
            pd.DataFrame.to_csv(dataToWrite, outputfilename, sep='\t', na_rep='NULL')
        else:
           with open(outputfilename, "wt") as out_file:
               out_file.write("No Significant Interactors")

    def GenerateKDEAlonePlots(self, outputFileName):
        ##Read in Data
        df = self._LFQDataFrame
        ##Read in Curve
        df2 = self._CurveDataFrame
        ## NOTES
        ## replace NaN in subcomplex labels with unlabeled so that it plots correctly, and then setup color
        ## hue so that it turns grey
        ## 
        df = df.replace({'SubComplexAnnotation': { np.NaN : 'Unlabeled'}})


        ## Think about doing a 1D density plot on the axes instead of a contour plot

        # Subset the Volcano Plot by ProteasomeSubComplex

        RP = df.loc[df['SubComplexAnnotation'] == 'RP']
        CP = df.loc[df['SubComplexAnnotation'] == 'CP']
        PIPS_CP = df.loc[df['SubComplexAnnotation'] == 'PIP_CP']
        PIPS_RP = df.loc[df['SubComplexAnnotation'] == 'PIP_RP']
        PIP = df.loc[df['SubComplexAnnotation'] == 'PIP']

        ## Setup XKCD Colors
        colors = ["cool grey", "windows blue", "pale red", "medium green", "orange", "black"]
        ## Parse Colors
        colors = sns.xkcd_palette(colors);

        ##SetupCusotmPalleteDictionary
        ##if self._LabelType=="AllPips":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[3], PIP_RP=colors[3])
        ##if self._LabelType=="CPPips":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[3], PIP_RP=colors[0])
        ##if self._LabelType=="RPPips":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[0], PIP_RP=colors[3])
        ##if self._LabelType=="ECM29":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[3], PIP_CP=colors[0], PIP_RP=colors[0])
        ##if self._LabelType=="CPRP":
        ##    custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[0], PIP_CP=colors[0], PIP_RP=colors[0])
    
        ##Rick's Palette
        custompalette=dict(Unlabeled=colors[0], RP=colors[1], CP=colors[2], PIP=colors[5], PIP_CP=colors[3], PIP_RP=colors[4])

        red = sns.color_palette("Reds")[-2]
        blue = sns.color_palette("Blues")[-2]

        ## Setup Seaborn facetgrid with lmplot (no line)
        plt.gcf().clear()
        kdeax = sns.kdeplot(CP['Difference'], label='CP Density', color =red, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5)
        kdemedianset1 = self.KDEMedian(kdeax)
        plt.gcf().clear()
        kdeax2 = sns.kdeplot(RP['Difference'], label='RP Density', color =blue, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5)
        kdemedianset2 = self.KDEMedian(kdeax2)
        plt.gcf().clear()

        kdeax = sns.kdeplot(CP['Difference'], label='CP Density', color =red, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5)
        kdeax2 = sns.kdeplot(RP['Difference'], label='RP Density', color =blue, bw=0.5, shade=True, shade_lowest=False, alpha = 0.5, ax =kdeax)

        ## Setup Legend
        #facetgrid.fig.get_axes()[0].legend(loc='upper left', title = 'Key')

        ## Setup Axes and axes limit
        myaxes = kdeax.axes
        myaxes.set_ylim(0,1)
        myaxes.set_xlim(-10,13)
        
        ## KS two sample test
        ks = scipy.stats.ks_2samp(CP['Difference'],RP['Difference'])
        mediantest = scipy.stats.median_test(CP['Difference'],RP['Difference'])
        mediantestpvalue = '%.2E' % Decimal(mediantest[1])
        mediantestpvalue = "Median P Value =" + str(mediantestpvalue)
        kspvalue = '%.2E' % Decimal(ks.pvalue)
        kspvalue = "KS P Value =" + str(kspvalue)


        plt.vlines(kdemedianset1[0], 0, kdemedianset1[1], color=red)
        plt.vlines(kdemedianset2[0], 0, kdemedianset2[1], color=blue)
        foldchange = np.abs(kdemedianset1[0] - kdemedianset2[0])
        foldchange = np.exp2(foldchange)
        foldchange = round(Decimal(foldchange), 2)
        xoffset = -5;
        yoffset = -0.2
        xposition = 13
        yposition = 1
        labeltoAdd = "Fold Change:"
        labeltoAdd += str(foldchange)
        plt.text(xposition+xoffset,yposition+yoffset, labeltoAdd, size=16)
        plt.text(xposition+xoffset,yposition+yoffset-0.05, kspvalue, size=16)
        plt.text(xposition+xoffset,yposition+yoffset-0.1, mediantestpvalue, size=16)


        ## Save Figure
        plt.savefig(outputFileName + ".pdf");
        plt.savefig(outputFileName + ".png");    
        

    def LabelFigure(self, dataframe, ax, offset, colors):
        for index, row in dataframe.iterrows():
            xposition=row['Difference']
            yposition=row['-Log(P-value)']
            labeltoAdd=row['ProteasomeAnnotation']
            ax.text(xposition+offset,yposition+offset, labeltoAdd, size=16, color = colors[5])
    
    def KDEMedian(self, kdeplot):
        x,y = kdeplot.get_lines()[0].get_data()
        cdf = scipy.integrate.cumtrapz(y, x, initial=0)
        nearest_05 = np.abs(cdf-0.5).argmin()

        x_median = x[nearest_05]
        y_median = y[nearest_05]

        return [x_median, y_median]
            
    def Thoughts(self):
    #LabelList = ["PAG1_Minus_vs_Col0_Minus",
    #    "PAG1_Minus_vs_PAG1_Plus",
    #    "PAG1_Plus_vs_Col0_Plus",
    #    "RPT4a_Minus_vs_Col0_Minus",
    #    "RPT4a_Minus_vs_RPT4b_Minus",
    #    "RPT4a_Plus_vs_Col0_Plus",
    #    "RPT4a_Plus_vs_RPT4b_Plus",
    #    "RPT4b_Minus_vs_Col0_Minus",
    #    "RPT4b_Minus_vs_RPT4b_Plus",
    #    "RPT4b_Plus_vs_Col0_Plus",
    #    "RPTa_Minus_vs_RPT4a_Plus"]


    #TypeList = ["CPPips",
    #            "CPPips",
    #            "CPPips",
    #            "RPPips",
    #            "RPPips",
    #            "RPPips",
    #            "RPPips",
    #            "RPPips",
    #            "RPPips"] Length here is wrong, make same dimensions as LabelList
        test = self

    def GenerateBarGraph(self, dataframe, outputFileName, GraphGroupName):
        dataToPlot = df.loc[dataframe['GraphGroup' == GraphGroupName]]



## General Script

fileAndStyle = {'PAG1_Minus_Col0_Minus': 'CPPips','PAG1_Plus_Col0_Plus': 'CPPips','RPT4a_Minus_Col0_Minus': 'RPPips', 'RPT4a_Plus_Col0_Plus': 'RPPips', 'RPT4b_Minus_Col0_Minus': 'RPPips', 'RPT4b_Plus_Col0_Plus': 'RPPips' }

fileAndStyle2 = {'RPT4b_Minus_RPT4a_Minus' : 'CPRP', 'RPT4b_Plus_RPT4a_Plus' : 'CPRP' }

fileAndStyle3 = {'RPT4b_Plus_RPT4b_Minus' : 'ECM29', 'RPT4a_Plus_RPT4a_Minus' : 'ECM29', 'PAG1_Plus_PAG1_Minus' : 'ECM29'}

for name, labelstyle in fileAndStyle.items():
    plot = PlottingFunctions(name +".txt", name + "_Curve.txt", labelstyle)
    plot.GenerateFigure(name, shouldLabel=True)
    plot.GenerateKDEAlonePlots(name + "_kde.txt")
    plot.GetSignificantInteractorsThatAreNotProteasomeSubunits(name + "_Significant_Interactors.txt")

for name, labelstyle in fileAndStyle2.items():
    plot = PlottingFunctions(name +".txt", name+ "_Curve.txt", labelstyle)
    plot.GenerateFigure(name, shouldLabel=True)
    plot.GenerateKDEAlonePlots(name + "_kde.txt")
    plot.GetSignificantInteractorsThatAreNotProteasomeSubunits(name + "_Significant_Interactors.txt")