import seaborn as sns
import matplotlib as mpl
import numpy as np
import DCG_Utilities as dcgutil
## Allows for easy text manipulation in illustrator by changing this
mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['svg.fonttype'] = 'none'

import matplotlib.pyplot as plt
import pandas as pd

font = {'family' : 'Arial',
        'weight' : 'bold',
        'size' : 14}
mpl.rc('font',**font)


## Setup figure Style
#sns.set(context="paper", font_scale=2, style = "white")
sns.set(context="paper", font_scale=2)
sns.despine()

dataFrameList = []

class BarGraphs():
    def __init__(self, **kwargs):
        return super().__init__(**kwargs)

    def GenerateBarGraph(self, dataframe, outputFileName, GraphGroupName):
        dataToPlot = dataframe.loc[dataframe['C: GraphGroup']==GraphGroupName]
        #dataToPlot = dataframe
        #dataToPlot = self.RestructureDataFrameForSeaborn(self, dataframe=dataToPlot)
        ColumnFilterSansCol0 =('PAG1_Plus','PAG1_Minus', 'RPT4a_Plus', 'RPT4a_Minus', 'RPT4b_Plus','RPT4b_Minus','C: ProteasomeAnnotation', 'C: GraphGroup', 'C: SubComplexAnnotation')
        ColumnFilter=('Col0_Plus','Col0_Minus','PAG1_Plus','PAG1_Minus', 'RPT4a_Plus', 'RPT4a_Minus', 'RPT4b_Plus','RPT4b_Minus','C: ProteasomeAnnotation', 'C: GraphGroup', 'C: SubComplexAnnotation')
        dataToPlot = dcgutil.Utilities.FilterDataFramebyInclusionList(dataframe=dataToPlot, inclusionlist =ColumnFilterSansCol0, sortbycolumnname ='C: ProteasomeAnnotation')
        dataToPlot = self.DataMassage(self, dataframe=dataToPlot)
        dataFrameList.append(dataToPlot)
        rpnlist = ["RPN1a", "RPN1b", "RPN2a", "RPN3a", "RPN3b", "RPN5a", "RPN5b",
                   "RPN6", "RPN7", "RPN8a", "RPN8b", "RPN9a", "RPN9b", "RPN10",
                   "RPN11", "RPN12a", "Sem1(DSS1)-1"]
        associatedlist = ["PBAC1", "PBAC2", "PBAC3", "PBAC4", "PAP1", "Ump1a",
                          "Ump1c?", "HSM3", "Nas2", "Nas6", "PA200", "ECM29", "PTRE1"]
        orderlookup ={'RPN' : rpnlist, 'Assembly' : associatedlist}
        self.PlotBarGraph(self,dataframe=dataToPlot, outputFileName=outputFileName, order=orderlookup.get(GraphGroupName))
        

    def DataMassage(self, dataframe):
        mydataframe = dataframe
        #print(mydataframe)
        #check = mydataframe['C: SubComplexAnnotation']
        #if 'Assembly' in check.values:
        #    mydataframe.drop['C: SubComplexAnnotation']
        mydataframe = mydataframe.set_index(['C: ProteasomeAnnotation', 'C: GraphGroup', 'C: SubComplexAnnotation'])
        mydataframe.index.names = ['Subunit', 'GraphGroup', 'SubComplex']
        ##Maybe need to index based on subunit
        #arrays = [np.array(['Col0', 'Col0', 'Col0', 'Col0', 'Col0', 'Col0',
        #                    'PAG1', 'PAG1', 'PAG1', 'PAG1', 'PAG1', 'PAG1',
        #                    'RPT4a','RPT4a','RPT4a','RPT4a','RPT4a','RPT4a',
        #                    'RPT4b','RPT4b','RPT4b','RPT4b','RPT4b','RPT4b',]),
        #          np.array(['Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',
        #                    'Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',
        #                    'Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',
        #                    'Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',])]

        ## For Col0 Plots
        #tuples = list(zip(*[['Col0', 'Col0', 'Col0', 'Col0', 'Col0', 'Col0',
        #                    'PAG1', 'PAG1', 'PAG1', 'PAG1', 'PAG1', 'PAG1',
        #                    'RPT4a','RPT4a','RPT4a','RPT4a','RPT4a','RPT4a',
        #                    'RPT4b','RPT4b','RPT4b','RPT4b','RPT4b','RPT4b',],
        #                    ['Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',
        #                    'Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',
        #                    'Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',
        #                    'Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',]]))
        ## For Just CP and RP Plots
        tuples = list(zip(*[['PAG1', 'PAG1', 'PAG1', 'PAG1', 'PAG1', 'PAG1',
                            'RPT4a','RPT4a','RPT4a','RPT4a','RPT4a','RPT4a',
                            'RPT4b','RPT4b','RPT4b','RPT4b','RPT4b','RPT4b',],
                            ['Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',
                            'Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',
                            'Minus', 'Minus', 'Minus', 'Plus', 'Plus', 'Plus',]]))


        newcolumns = pd.MultiIndex.from_tuples(tuples, names= ['Pulldown','ATP'])
        mydataframe.columns = newcolumns
        mydataframe = mydataframe.T
        columnnames = list(mydataframe.index.levels[0])
        mydataframe = mydataframe.stack('Subunit')
        print(mydataframe)
        mydataframe = mydataframe.stack('SubComplex')
        mydataframe = mydataframe.stack('GraphGroup')
        print(mydataframe)
        mydataframe = mydataframe.reset_index(level=['Pulldown','ATP','Subunit', 'SubComplex', 'GraphGroup'])
        print(mydataframe)
        mydataframe.columns= ['Pulldown','ATP','Subunit','SubComplex','GraphGroup','Log2(LFQ)']
        print(mydataframe)
        mydataframe['PulldownLabel'] = mydataframe['Pulldown'] + mydataframe['ATP']
        return mydataframe
    
    def PlotBarGraph(self, dataframe, outputFileName, order = None):
        ## Setup XKCD Colors
        
        #labelorderdict = { "Alpha" : "test"};

        colors = ["cool grey", "pale red", "windows blue", "light blue"]
        #rickscolors = ["#000000","#fbb016","#7270b3","#b3d234"]
        rickscolors = [(255/255,255/255,255/255),(251/255,176/255,22/255),(114/255,112/255,179/255),(179/255,210/255,52/255)]

        ## Parse Colors
        #rickscolors = sns.color_palette(rickscolors)
        #colors = sns.xkcd_palette(colors);
        ##SetupCusotmPalleteDictionary
        custompalette=dict(Col0Minus=rickscolors[0], rickscolors=colors[0],
                           PAG1Minus=rickscolors[1], PAG1Plus=rickscolors[1],
                           RPT4aMinus=rickscolors[2], RPT4aPlus=rickscolors[2],
                           RPT4bMinus=rickscolors[3], RPT4bPlus=rickscolors[3])

        if(order != None):
            fg = sns.factorplot(x='Subunit', y='Log2(LFQ)', hue='PulldownLabel', kind='bar', col='ATP',
                           data=dataframe, size=12, palette=custompalette, conf_lw=1.5, capsize=0.1,
                           ci=68, order=order, saturation=1)
        else:
            fg = sns.factorplot(x='Subunit', y='Log2(LFQ)', hue='PulldownLabel', kind='bar', col='ATP',
                           data=dataframe, size=12, palette=custompalette, conf_lw=1.5, capsize=0.1,
                           ci=68, saturation = 1)

        #fg = sns.factorplot(x='Subunit', y='Log2(LFQ)', hue='PulldownLabel', kind='bar', col='ATP', row = 'SubComplex',
        #                    data=dataframe, size=12, palette=custompalette, conf_lw=1.5, capsize=0.1, ci=68)
        fg.set(ylim=(15,35))
        fg.set_xticklabels(rotation=-30)
        fg.savefig(outputFileName + ".png")
        fg.savefig(outputFileName + ".pdf")

    def PlotBarGraphs(self, dataframeList, outputFileName):
        ## Setup XKCD Colors
        colors = ["cool grey", "pale red", "windows blue", "light blue"]
        ## Parse Colors
        colors = sns.xkcd_palette(colors);
        ##SetupCusotmPalleteDictionary
        custompalette=dict(Col0Minus=colors[0], Col0Plus=colors[0],
                           PAG1Minus=colors[1], PAG1Plus=colors[1],
                           RPT4aMinus=colors[2], RPT4aPlus=colors[2],
                           RPT4bMinus=colors[3], RPT4bPlus=colors[3])

        ##ConcatenateDataFrame
        concatDataFrame = pd.DataFrame()
        for dataframes in dataframeList:
            concatDataFrame = pd.concat([concatDataFrame, dataframes])
        print(concatDataFrame)

        #Depreciated
        #concatDataFrame = concatDataFrame.assign(GraphGroup=concatDataFrame.Subunit.astype(object)).sort("Subunit")
        #hueorder = concatDataFrame.PulldownLabel.unique()
        #fg = fg.map_dataframe(sns.barplot,"Subunit", "Log2(LFQ)", hue="PulldownLabel")
        #fg = sns.factorplot('Subunit', y='Log2(LFQ)', hue='PulldownLabel', kind='bar', col='ATP', row = 'GraphGroup',
        #                    data=concatDataFrame, size=12, palette=custompalette, conf_lw=1.5, capsize=0.1, ci=68, sharex = 'row')
        
        ##Setup FacetGrid with Plus and Minus ATP and Row Equal to GraphGroup
        fg = sns.FacetGrid(data = concatDataFrame, col="ATP", row="GraphGroup",
                           sharex = False, size =12)
        ##Map custom function MyBarPlot to give good defaults onto the facetgrid
        fg = fg.map_dataframe(self.MyBarPlot, "Subunit", "Log2(LFQ)","PulldownLabel",
                              custompalette)
        ##Get All x axes in the Facetgrid and set their rotation to -30
        for i, ax in enumerate(fg.axes.flat):
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=-30)
        ##Set ylim to 15
        fg.set(ylim=(15, None))

        #SaveFgiures
        fg.savefig(outputFileName + ".png")
        fg.savefig(outputFileName + ".pdf")

    def MyBarPlot(x, y, hue=None, palette=None, data=None, label=None, color=None, **kwargs):
        sns.barplot(data = data, x=x,y=y,hue=hue, palette=palette, conf_lw=1.5, capsize=0.1, ci=68)
    def StandardizeColumns(self, dataframe):
        dataframe.reset_index(inplace=True)
        for column in dataframe:
            if "RPT4a" in column:
                RPT4avalue = np.where(dataframe['C: ProteasomeAnnotation']=="RPT4a")[0].tolist()[0]
                RPT4avalue = dataframe.ix[RPT4avalue, column]
                dataframe[column]= dataframe[column].apply(lambda x : x / RPT4avalue)
            if "RPT4b" in column:
                RPT4bvalue = np.where(dataframe['C: ProteasomeAnnotation']=="RPT4b")[0].tolist()[0]
                RPT4bvalue = dataframe.ix[RPT4bvalue, column]
                dataframe[column]= dataframe[column].apply(lambda x : x / RPT4bvalue)
            if "PAG1" in column:
                PAG1value = np.where(dataframe['C: ProteasomeAnnotation']=="PAG1")[0].tolist()[0]
                PAG1value = dataframe.ix[PAG1value, column]
                dataframe[column]= dataframe[column].apply(lambda x : x / PAG1value)
        return dataframe;
    def StandardizeColumnsRPT4aRPT4b(self, dataframe):
        dataframe.reset_index(inplace=True)
        #for column in dataframe:
        #    if "RPT4a" in column or "RPT4b" in column:
        #        RPT4avalue = np.where(dataframe['C: ProteasomeAnnotation']=="RPT4a")[0].tolist()[0]
        #        RPT4avalue = dataframe.ix[RPT4avalue, column]
        #        RPT4bvalue = np.where(dataframe['C: ProteasomeAnnotation']=="RPT4b")[0].tolist()[0]
        #        RPT4bvalue = dataframe.ix[RPT4bvalue, column]
        #        RPT4a_corrected = (RPT4avalue * (RPT4bvalue
        #        RPT4b_corrected =


bargraph = BarGraphs
GraphList = ['Alpha','Beta','RPT','RPN', 'Assembly']
#GraphList = ['Assembly']

for graphgroup in GraphList:
    bargraph.GenerateBarGraph(bargraph,
                              dataframe = pd.read_table("Bar_Graph_No_Imputation_All_BioReps.txt").fillna(value = 15),
                              #dataframe = bargraph.StandardizeColumns(bargraph, pd.read_table("Linear_Scale_Data_Bar_Graphs.txt")),
                              outputFileName= "BarGraph" + graphgroup, GraphGroupName= graphgroup)

#bargraph.PlotBarGraphs(bargraph, dataframeList=dataFrameList, outputFileName="All")