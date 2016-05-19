import seaborn as sns
import matplotlib as mpl
import numpy as np
from decimal import *
from scipy.stats import linregress
## Allows for easy text manipulation in illustrator by changing this
mpl.rcParams['pdf.fonttype'] = 42

import matplotlib.pyplot as plt
import pandas as pd
import DCG_Utilities as dcgutil

font = {'family' : 'Arial',
        'weight' : 'bold',
        'size' : 14}
mpl.rc('font',**font)


## Setup figure Style
sns.set(context="paper", font_scale=2, style = "white")
#sns.set(context="paper", font_scale=2)
sns.despine()

dataFrameList = []

class ScatterGraphs:
    def __init__(self, **kwargs):
        return super().__init__(**kwargs)

    def ProcessData(self, dataframe):
        dataframeToReturn = pd.DataFrame();
        ColumnFilter=('Col0_Plus','Col0_Minus','PAG1_Plus','PAG1_Minus', 'RPT4a_Plus', 'RPT4a_Minus', 'RPT4b_Plus','RPT4b_Minus','C: ProteasomeAnnotation')
        dataframeToReturn = dcgutil.Utilities.FilterDataFramebyInclusionList(dataframe,ColumnFilter, "C: ProteasomeAnnotation")
        
        del dataframeToReturn['C: ProteasomeAnnotation']
        columnames=list(dataframeToReturn.columns)
        dataframeToReturn = dataframeToReturn.unstack()
        dataframeToReturn = pd.Series.to_frame(dataframeToReturn)
        dataframeToReturn.reset_index(level=0, inplace=True)
        dataframeToReturn.columns = ['Label', 'LFQ']
        dataframeToReturn['Pulldown'] = None
        dataframeToReturn['ATP'] = None
        dataframeToReturn['PulldownLabel'] = None
        dataframeToReturn['BioRep'] = None
        dataframeToReturn['TechRep'] = None
        dataframeToReturn = dataframeToReturn.reset_index()
        print(dataframeToReturn)
        for i in dataframeToReturn.index:
            string = dataframeToReturn.iloc[i]['Label']
            #string = series.iloc[i]
            strings = str.split(string)
            stringsToSplitAgain = strings[len(strings)-1]
            strings = str.split(stringsToSplitAgain,'_')
            dataframeToReturn.set_value(i,'Pulldown', strings[0])
            dataframeToReturn.set_value(i,'ATP', strings[1])
            dataframeToReturn.set_value(i,'PulldownLabel', strings[0] + "_" + strings[1])
            dataframeToReturn.set_value(i,'BioRep', strings[2])
            dataframeToReturn.set_value(i,'TechRep', strings[3])
            print(i)

        print(dataframeToReturn)
        #for index in range(0,len(columnames)-1,2):
        #    data1 = columnames[index]
        #    data2 = columnames[index+1]
        #    figure = sns.lmplot(data1, data2, data=dataframeToReturn)
        #    figure.savefig(columnames[index] + columnames[index+1] + ".png")
            
        dataframeToReturn.to_pickle("temp.pickle")
        return dataframeToReturn
    
    def PlotScatterPlots(self, dataframe):
        colors = ["cool grey", "pale red", "windows blue", "light blue"]
        colors = sns.xkcd_palette(colors)
        #dataframe["LFQ_TR1"] = dataframe.loc[dataframe["TechRep"] =="TR1", "LFQ"]
        #dataframe["LFQ_TR2"] = dataframe.loc[dataframe["TechRep"] =="TR2", "LFQ"]
        #print(dataframe)
        #print(dataframe)

        fg = sns.FacetGrid(data = dataframe, col="BioRep", row="PulldownLabel", size =12, sharex=False, sharey=False)

        fg = fg.map_dataframe(self.MyScatterPlot, "LFQ", "LFQ")

        #sns.FacetGrid.map_dataframe

        fg.savefig("scatterTest.png")
        fg.savefig("scatterTest.pdf")

    def MyScatterPlot(x, y, label=None, color=None, **kwargs):
        
        colors = ["cool grey", "pale red", "windows blue", "light blue"]
        colors = sns.xkcd_palette(colors)

        data = kwargs.pop('data')
        data1 = data.loc[data["TechRep"] =="TR1"]
        data2 = data.loc[data["TechRep"] =="TR2"]

        #How to select subset of data
        #data3 = data.loc[data["TechRep"] =="TR1", "LFQ"]

        x = data1["LFQ"]
        y = data2["LFQ"]
        #data1 = data.loc[data["TechRep"] == "TR1"]
        #data1 = data.loc[data["LFQ"]]
        #data2 = data.loc[data["TechRep"] == "TR2"]
        #data2 = data.loc[data["LFQ"]]
        #print(data1)
        #print(data2)
        

        ##Calculate Regression Coeficients
        reg = linregress(x,y)
        slope = reg.slope
        intercept = reg.intercept
        rvalue = reg.rvalue
        rsquared = reg.rvalue * reg.rvalue
        rsquared = round(Decimal(rsquared), 3)
        #sns.regplot(

        myfigure = plt.scatter(x=x, y=y, s = 200, color="grey")
        
        


        ## Regresion
        plt.plot(x, slope*x + intercept, color = "black")
        ## Add R^2 Text
        plt.annotate("R^2=" + str(rsquared), (0.5, 0.9), textcoords='axes fraction', size=40)
        plt.xlim(xmin=-0.1)
        plt.ylim(ymin=-0.1)
        

        #figure = sns.lmplot(x, y, kwargs.pop('data'), palette = colors)
        #figure.savefig("wtf.pdf")

        print("figure")

    def ExperimentalRearrange(dataframe):
        
        print(dataframe)
        dataframe = dataframe.reset_index()
        dataframe = dataframe.pivot(columns= "TechRep", values = "LFQ")

        print(dataframe)

        dataframe = pd.melt(dataframe, id_vars ="TechRep", value_vars ="LFQ")

        print(dataframe)
            
        return dataframe
        

scattergraph = ScatterGraphs

scattergraph.ProcessData(scattergraph, pd.read_table("Tech_Rep_Correlations.txt"))

scattergraph.PlotScatterPlots(scattergraph, pd.read_pickle("temp.pickle"))

#scattergraph.ExperimentalRearrange(pd.read_pickle("temp.pickle"))
#print(mydataframe)


##Row should be Sample Label
##Col should be bio Rep
##X should be tech rep 1, Y should be tech rep 2