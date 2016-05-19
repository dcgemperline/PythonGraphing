import seaborn as sns
import matplotlib as mpl
import numpy as np
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


class BarGraphs():
    def __init__(self, **kwargs):
        return super().__init__(**kwargs)

    def GenerateBarGraph(self, dataframe, outputFileName, GraphGroupName):
        dataToPlot = dataframe.loc[dataframe['C: GraphGroup']==GraphGroupName]
        #dataToPlot = dataframe
        #dataToPlot = self.RestructureDataFrameForSeaborn(self, dataframe=dataToPlot)
        ColumnFilterSansCol0 =('PAG1_Plus','PAG1_Minus', 'RPT4a_Plus', 'RPT4a_Minus', 'RPT4b_Plus','RPT4b_Minus','C: ProteasomeAnnotation', 'C: GraphGroup', 'C: SubComplexAnnotation')
        ColumnFilter=('Col0_Plus','Col0_Minus','PAG1_Plus','PAG1_Minus', 'RPT4a_Plus', 'RPT4a_Minus', 'RPT4b_Plus','RPT4b_Minus','C: ProteasomeAnnotation', 'C: GraphGroup', 'C: SubComplexAnnotation')
        dataToPlot = self.RemoveUnneccesaryColumns(self, dataframe=dataToPlot, inclusionlist =ColumnFilterSansCol0, sortbycolumnname ='C: ProteasomeAnnotation')
        dataToPlot = self.DataMassage(self, dataframe=dataToPlot)
        self.PlotBarGraph(self,dataframe=dataToPlot, outputFileName=outputFileName)


    def RemoveUnneccesaryColumns(self, dataframe, inclusionlist, sortbycolumnname):
        df = dataframe
        columnnames = list(df)
        for names in columnnames:
            found = False
            for inclusionname in inclusionlist:
                if inclusionname in names:
                    found = True
            if not found == True:
                del df[names]
        df = df.sort_values(by=sortbycolumnname, ascending=True)
        return df

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
        ## For Just CP Plots
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
        mydataframe = mydataframe.stack('SubComplex')
        mydataframe = mydataframe.reset_index(level=['Pulldown','ATP','Subunit', 'SubComplex'])
       # if 'Assembly' in columnnames:
            #mydataframe = mydataframe.unstack('SubComplex')
            #mydataframe = mydataframe.dropna()
        print(mydataframe)
        mydataframe.columns= ['Pulldown','ATP','Subunit','SubComplex','Log2(LFQ)']
        print(mydataframe)
        mydataframe['PulldownLabel'] = mydataframe['Pulldown'] + mydataframe['ATP']
        return mydataframe

    def GenerateLists(stringArray1, stringArray2, label1num, label2repeatnum):
        ## Generate Tuples for the zip function to make dynamic analysis for specific things easier
        
        mylist1 = []
        mylist2 = []
        for string in range(0, len(stringArray1)-1):
            for iter in range(0,label1num-1):
                mylist1.append(string)
        for iter2 in range(0, (len(mylist1)-1)/ label2repeatnum):
            for iter2 in range(0, label2repeatnum-1):
                test = test
        


    def PlotBarGraph(self, dataframe, outputFileName):
        ## Setup XKCD Colors
        colors = ["cool grey", "pale red", "windows blue", "light blue"]
        ## Parse Colors
        colors = sns.xkcd_palette(colors);
        dataframe = dataframe
        ##SetupCusotmPalleteDictionary
        custompalette=dict(Col0Minus=colors[0], Col0Plus=colors[0],
                           PAG1Minus=colors[1], PAG1Plus=colors[1],
                           RPT4aMinus=colors[2], RPT4aPlus=colors[2],
                           RPT4bMinus=colors[3], RPT4bPlus=colors[3])

        #fg = sns.factorplot(x='Subunit', y='Log2(LFQ)', hue='PulldownLabel', kind='bar', col='ATP',
        #                   data=dataframe, size=12, palette=custompalette, conf_lw=1.5, capsize=0.1, ci=68)

        fg = sns.factorplot(x='Subunit', y='Log2(LFQ)', hue='PulldownLabel', kind='bar', col='ATP', row = 'SubComplex',
                            data=dataframe, size=12, palette=custompalette, conf_lw=1.5, capsize=0.1, ci=68)
        fg.set(ylim=(15, None))
        fg.set_xticklabels(rotation=-30)
        fg.savefig(outputFileName + ".png")
        fg.savefig(outputFileName + ".pdf")

           
    def RestructureDataFrameForSeaborn(self, dataframe):
        ## Don't use this it doesnt work
        
        df = pd.DataFrame()
        tempdf = pd.DataFrame()
        dftoreturn = pd.DataFrame()
        
        df = dataframe
        ##SubsetDataFrame

        ColumnFilter=('Col0_Plus','Col0_Minus','PAG1_Plus','PAG1_Minus', 'RPT4a_Plus', 'RPT4a_Minus', 'RPT4b_Plus','RPT4b_Minus','C: ProteasomeAnnotation', 'C: GraphGroup','C: ProteasomeAnnotation', 'C: SubComplexAnnotation')
        first = True
        for string in ColumnFilter:
            tempdf = df.filter(regex=string)
            if first:
                dftoreturn = tempdf
                first = False
            else:
                for columns in tempdf:
                    dftoreturn.append(tempdf[columns])
        dftoreturn = dftoreturn.sort_values(by='C: ProteasomeAnnotation', ascending=True)
        return dftoreturn


bargraph = BarGraphs
GraphList = ['Alpha','Beta','RPT','RPN', 'Assembly']
#GraphList = ['Assembly']

for graphgroup in GraphList:
        bargraph.GenerateBarGraph(bargraph,
                          dataframe = pd.read_table("Bar_Graph_Data_Imputed_All_BioReps.txt"),
                          outputFileName= "BarGraph" + graphgroup, GraphGroupName= graphgroup)

