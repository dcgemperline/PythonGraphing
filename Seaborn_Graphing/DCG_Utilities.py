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


class Utilities:

    @staticmethod
    def FilterDataFramebyInclusionList(dataframe, inclusionlist, sortbycolumnname):
        columnnames = list(dataframe)
        for names in columnnames:
            found = False
            for inclusionname in inclusionlist:
                if inclusionname in names:
                    found = True
            if not found == True:
                del dataframe[names]
        dataframe = dataframe.sort_values(by=sortbycolumnname, ascending=True)
        return dataframe
