import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns

#df1 = pd.read_table('raw_FlashQE.txt', skiprows=11) #skipping first 11 rows
#EQE=df1[(df1['Param'])=='EQE']

in_data=pd.read_table('data.txt', skiprows=0)

'''
#cols = ['Isc_(A)', 'Voc_(V)', 'FF_(percent)', 'pJmp_(A/cm2)', 'pVmp_(V)', 'pPmp_(W/cm2)', 'n_at_1_sun', 'n_at_ 1/10_suns', 'Jo1_(A/cm2)', 'Jo2_(A/cm2)', 'Lifetime_at_Vmp_(us)', 'Doping_(cm-3)', 'Measured_Resistivity_(Ohm-cm)', 'Capacitance_(uF/cm2)', 'BulkSRH_AND_Jo_Power_Loss_(percent)', 'Rsh_Power_Loss_(percent)', 'AverageJsc', 'DeadLayerThickness_Average', 'EffectiveDiffusionLength_Average', 'ReflectionMinWavelength_Average', 'FrontSurfaceReflectanceLoss_Average', 'EscapeReflectanceLoss_Average', 'EmitterLoss_Average', 'BulkandRearLoss_Average', 'Mean V_xy']
#sns.pairplot(EQE[cols], size=2.5)
cols = ['Isc_(A)', 'Voc_(V)', 'FF_(percent)']
g = sns.FacetGrid(EQE[cols])
'''

g = sns.PairGrid(in_data)
g = g.map(plt.scatter)

xlabels,ylabels = [],[]

for ax in g.axes[-1,:]:
    xlabel = ax.xaxis.get_label_text()
    xlabels.append(xlabel)
for ax in g.axes[:,0]:
    ylabel = ax.yaxis.get_label_text()
    ylabels.append(ylabel)

for i in range(len(xlabels)):
    for j in range(len(ylabels)):
        g.axes[j,i].xaxis.set_label_text(xlabels[i])
        g.axes[j,i].yaxis.set_label_text(ylabels[j])

plt.show()
g.savefig('pvsc_big_data.jpg', format='jpg', dpi=400)
