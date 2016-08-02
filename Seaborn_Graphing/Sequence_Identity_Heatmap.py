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
sns.set(context="paper", font_scale=1, style = "white")
sns.despine()

proteins = ["pap1", "pbac1", "pbac2", "pbac3", "pbac4", "ump1"]

for protein in proteins:
    plt.gcf().clear()
    df = pd.read_table(protein + ".ident.matrix", skiprows=1, header=None, delim_whitespace=True, index_col=0)
    names = df.index.tolist()
    df.index.name = "Sequences"
    df.columns = names
    mask = np.zeros_like(df)
    mask[np.triu_indices_from(mask)] = True
    ax = sns.heatmap(data=df, square=True, vmin=0, vmax=100, mask=mask, annot=True, fmt='.1f', cmap="Blues")
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    ax.figure.savefig(protein + "identity_matrix.pdf")
    