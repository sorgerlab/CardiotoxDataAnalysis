#Scatter plots
from sklearn.cluster import KMeans
from statistics import mode
from statistics import median
from statistics import mean
from statistics import StatisticsError
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Manual gates set from looking at images in Columbus
mlc2a_gate = 11.2
pRb_gate = 10.4
mlc2v_gate = 10.9
aact_gate = 10.9    #From clustering
#plt.axvline(x_position)



for ab in ['proteome1']:
    for plate in ['P6','P7','P8']:

#        for dose in ['lowdose','highdose']:
        df1 = pd.read_csv('%s/dmso/%s_dmso_ProcessedData%s_%s.csv' % ('lowdose','lowdose',plate,ab))
        df2 = pd.read_csv('%s/dmso/%s_dmso_ProcessedData%s_%s.csv' % ('highdose','highdose',plate,ab))
        df = pd.concat([df1,df2])


        for obs in ['Mlc2a_Cyto','MLC2v_Cyto','a-ACTININ_Cyto','mlc_total']:
            mlc2aHi = []
            mlc2aLo = []

            if obs == 'Mlc2a_Cyto':
                gate = mlc2a_gate
                MLC2aHiCells = df.loc[df[obs] >= gate]
                MLC2aLoCells = df.loc[df[obs] < gate]
                plt=df.plot.scatter('p-Rb_Nuc',obs)
                plt.axvline(pRb_gate)
                plt.axhline(gate)
                fig = plt.figure
                fig.savefig('%s_%s_p-Rb.png' % (plate,obs.split('_')[0]))


            elif obs == 'MLC2v_Cyto':
                gate = mlc2v_gate
                MLC2aHiCells = df.loc[df[obs] >= gate]
                MLC2aLoCells = df.loc[df[obs] < gate]
                plt=df.plot.scatter('p-Rb_Nuc',obs)
                plt.axvline(pRb_gate)
                plt.axhline(gate)
                fig = plt.figure
                fig.savefig('%s_%s_p-Rb.png' % (plate,obs.split('_')[0]))

            elif obs == 'ACTININ_Cyto':
                gate = aact_gate
                MLC2aHiCells = df.loc[df[obs] >= gate]
                MLC2aLoCells = df.loc[df[obs] < gate]
                plt=df.plot.scatter('p-Rb_Nuc',obs)
                plt.axvline(pRb_gate)
                plt.axhline(gate)
                fig = plt.figure
                fig.savefig('%s_%s_p-Rb.png' % (plate,obs.split('_')[0]))

            elif obs == 'mlc_total':
                MLC2aHiCells = df.loc[((df['Mlc2a_Cyto'] >= mlc2a_gate) | (df['MLC2v_Cyto'] >= mlc2v_gate))]
                MLC2aLoCells = df.loc[((df['Mlc2a_Cyto'] < mlc2a_gate) & (df['MLC2v_Cyto'] < mlc2v_gate))]
                mlc_max = df[['Mlc2a_Cyto','MLC2v_Cyto']].max(axis=1)
                df['mlc_max'] = mlc_max
                plt=df.plot.scatter('p-Rb_Nuc','mlc_max')
                plt.axvline(pRb_gate)
                plt.axhline(mlc2v_gate)
                fig = plt.figure
                fig.savefig('%s_%s_p-Rb.png' % (plate,obs))



            pRbLo_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] >= pRb_gate]

            pRbLo_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] >= pRb_gate]

            mlc2aHi.append(pRbLo_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aHi.append(pRbHi_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            mlc2aLo.append(pRbLo_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aLo.append(pRbHi_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            print(plate)
            print(obs)
            print(mlc2aHi)
            print(mlc2aLo)





