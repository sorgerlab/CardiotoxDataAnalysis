#Scatter plots
from sklearn.cluster import KMeans
from statistics import mode
from statistics import median
from statistics import mean
from statistics import StatisticsError
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def gate_on_var(series):
    percentHi_all = []

    for i in range(0,500):
        model = KMeans(n_clusters=2,init='k-means++')
#        model.fit(df[var].values.reshape(-1,1))
        model.fit(series[0].values.reshape(-1,1))
        labels = model.labels_  
        clusters = pd.DataFrame(data=labels,columns=['cluster'])
        centers = model.cluster_centers_
        #print(centers)

        #Identify hi and lo clusters
        hi_label = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
        lo_label = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

        #Split into var_hi df
        idxHi = clusters.index[clusters['cluster']==hi_label].tolist()
        hiCells = series.loc[idxHi]
        hiCells = hiCells.dropna()
        idxLo = clusters.index[clusters['cluster']==lo_label].tolist()
        loCells = series.loc[idxLo]
        loCells = loCells.dropna()

    threshold = mean([min(hiCells[0]),max(loCells[0])])

    return threshold


#Manual gates set from looking at images in Columbus
mlc2a_gate = 11.2
pRb_gate = 10.4
oct4a_gate = 12.5
ki76_gate = 9.2285068600382314
mlc2v_gate = 10.9


for ab in ['proteome1']:
    for plate in ['P7']:

#        for dose in ['lowdose','highdose']:
        for drug in ['dmso','soraf']:

            df1 = pd.read_csv('%s/%s/%s_%s_ProcessedData%s_%s.csv' % ('lowdose',drug,'lowdose',drug,plate,ab))
            df2 = pd.read_csv('%s/%s/%s_%s_ProcessedData%s_%s.csv' % ('highdose',drug,'highdose',drug,plate,ab))
            df = pd.concat([df1,df2])

            df.loc[df['p-Rb_Nuc'] == max(df['p-Rb_Nuc'])] = np.NaN
            df.loc[df['Ki67_Nuc'] == max(df['Ki67_Nuc'])] = np.NaN
            df.loc[df['MLC2v_Cyto'] == max(df['MLC2v_Cyto'])] = np.NaN
            df.loc[df['Mlc2a_Cyto'] == max(df['Mlc2a_Cyto'])] = np.NaN
            df.loc[df['Oct4a_Nuc'] == max(df['Oct4a_Nuc'])] = np.NaN

#            df3 = df.loc[df['p-Rb_Nuc'] == max(df['p-Rb_Nuc'])]


            for obs in ['Oct4a_Nuc','mlc_total']:

                if obs == 'Oct4a_Nuc':
                    mlcHi = []
                    mlcLo = []
                    oct4aHi = []
                    oct4aLo = []
                    gate = oct4a_gate
                    Oct4aHiCells = df.loc[df[obs] >= gate]
                    Oct4aLoCells = df.loc[df[obs] < gate]
                    plt=df.plot.scatter('p-Rb_Nuc',obs)
                    plt.axvline(pRb_gate)
                    plt.axhline(gate)
                    plt.set_title(drug.upper())
                    fig = plt.figure
                    fig.savefig('%s_%s_p-Rb_%s.png' % (plate,obs.split('_')[0],drug),dpi=600)

                    pRbLo_Oct4aHi_Cells = Oct4aHiCells.loc[Oct4aHiCells['p-Rb_Nuc'] < pRb_gate]
                    pRbHi_Oct4aHi_Cells = Oct4aHiCells.loc[Oct4aHiCells['p-Rb_Nuc'] >= pRb_gate]

                    pRbLo_Oct4Lo_Cells = Oct4aLoCells.loc[Oct4aLoCells['p-Rb_Nuc'] < pRb_gate]
                    pRbHi_Oct4Lo_Cells = Oct4aLoCells.loc[Oct4aLoCells['p-Rb_Nuc'] >= pRb_gate]

                    oct4aHi.append(pRbLo_Oct4aHi_Cells.shape[0]/(pRbLo_Oct4aHi_Cells.shape[0] + pRbHi_Oct4aHi_Cells.shape[0] + pRbLo_Oct4Lo_Cells.shape[0] + pRbHi_Oct4Lo_Cells.shape[0]))
                    oct4aHi.append(pRbHi_Oct4aHi_Cells.shape[0]/(pRbLo_Oct4aHi_Cells.shape[0] + pRbHi_Oct4aHi_Cells.shape[0] + pRbLo_Oct4Lo_Cells.shape[0] + pRbHi_Oct4Lo_Cells.shape[0]))

                    oct4aLo.append(pRbLo_Oct4Lo_Cells.shape[0]/(pRbLo_Oct4aHi_Cells.shape[0] + pRbHi_Oct4aHi_Cells.shape[0] + pRbLo_Oct4Lo_Cells.shape[0] + pRbHi_Oct4Lo_Cells.shape[0]))
                    oct4aLo.append(pRbHi_Oct4Lo_Cells.shape[0]/(pRbLo_Oct4aHi_Cells.shape[0] + pRbHi_Oct4aHi_Cells.shape[0] + pRbLo_Oct4Lo_Cells.shape[0] + pRbHi_Oct4Lo_Cells.shape[0]))

                    print(plate)
                    print(drug)
                    print(obs)
                    print(oct4aHi)
                    print(oct4aLo)



                    gate = oct4a_gate
                    Oct4aHiCells = df.loc[df[obs] >= gate]
                    Oct4aLoCells = df.loc[df[obs] < gate]
                    plt=df.plot.scatter('Ki67_Nuc',obs)
                    plt.axvline(ki76_gate)
                    plt.axhline(gate)
                    plt.set_title(drug.upper())
                    fig = plt.figure
                    fig.savefig('%s_%s_Ki67_%s.png' % (plate,obs.split('_')[0],drug),dpi=600)

                    mlcHi = []
                    mlcLo = []
                    oct4aHi = []
                    oct4aLo = []

                    ki67Lo_Oct4aHi_Cells = Oct4aHiCells.loc[Oct4aHiCells['Ki67_Nuc'] < ki76_gate]
                    ki67Hi_Oct4aHi_Cells = Oct4aHiCells.loc[Oct4aHiCells['Ki67_Nuc'] >= ki76_gate]

                    ki67Lo_Oct4aLo_Cells = Oct4aLoCells.loc[Oct4aLoCells['Ki67_Nuc'] < ki76_gate]
                    ki67Hi_Oct4aLo_Cells = Oct4aLoCells.loc[Oct4aLoCells['Ki67_Nuc'] >= ki76_gate]

                    oct4aHi.append(ki67Lo_Oct4aHi_Cells.shape[0]/(ki67Lo_Oct4aHi_Cells.shape[0] + ki67Hi_Oct4aHi_Cells.shape[0] + ki67Lo_Oct4aLo_Cells.shape[0] + ki67Hi_Oct4aLo_Cells.shape[0]))
                    oct4aHi.append(ki67Hi_Oct4aHi_Cells.shape[0]/(ki67Lo_Oct4aHi_Cells.shape[0] + ki67Hi_Oct4aHi_Cells.shape[0] + ki67Lo_Oct4aLo_Cells.shape[0] + ki67Hi_Oct4aLo_Cells.shape[0]))

                    oct4aLo.append(ki67Lo_Oct4aLo_Cells.shape[0]/(ki67Lo_Oct4aHi_Cells.shape[0] + ki67Hi_Oct4aHi_Cells.shape[0] + ki67Lo_Oct4aLo_Cells.shape[0] + ki67Hi_Oct4aLo_Cells.shape[0]))
                    oct4aLo.append(ki67Hi_Oct4aLo_Cells.shape[0]/(ki67Lo_Oct4aHi_Cells.shape[0] + ki67Hi_Oct4aHi_Cells.shape[0] + ki67Lo_Oct4aLo_Cells.shape[0] + ki67Hi_Oct4aLo_Cells.shape[0]))

                    print(plate)
                    print(drug)
                    print(obs)
                    print(oct4aHi)
                    print(oct4aLo)


                elif obs == 'mlc_total':
                    MLC2aHiCells = df.loc[((df['Mlc2a_Cyto'] >= mlc2a_gate) | (df['MLC2v_Cyto'] >= mlc2v_gate))]
                    MLC2aLoCells = df.loc[((df['Mlc2a_Cyto'] < mlc2a_gate) & (df['MLC2v_Cyto'] < mlc2v_gate))]
                    mlc_max = df[['Mlc2a_Cyto','MLC2v_Cyto']].max(axis=1)
                    df['mlc_max'] = mlc_max
                    df.loc[df['mlc_max'] == max(df['mlc_max'])] = np.NaN
                    plt=df.plot.scatter('p-Rb_Nuc','mlc_max')
                    plt.axvline(pRb_gate)
                    plt.axhline(mlc2v_gate)
                    plt.set_title(drug.upper())
                    fig = plt.figure
                    fig.savefig('%s_%s_p-Rb_%s.png' % (plate,obs,drug),dpi=600)

                    mlcHi = []
                    mlcLo = []
                    oct4aHi = []
                    oct4aLo = []

                    pRbLo_mlcHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] < pRb_gate]
                    pRbHi_mlcHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] >= pRb_gate]

                    pRbLo_mlcLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] < pRb_gate]
                    pRbHi_mlcLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] >= pRb_gate]

                    mlcHi.append(pRbLo_mlcHi_Cells.shape[0]/(pRbLo_mlcHi_Cells.shape[0] + pRbHi_mlcHi_Cells.shape[0] + pRbLo_mlcLo_Cells.shape[0] + pRbHi_mlcLo_Cells.shape[0]))
                    mlcHi.append(pRbHi_mlcHi_Cells.shape[0]/(pRbLo_mlcHi_Cells.shape[0] + pRbHi_mlcHi_Cells.shape[0] + pRbLo_mlcLo_Cells.shape[0] + pRbHi_mlcLo_Cells.shape[0]))

                    mlcLo.append(pRbLo_mlcLo_Cells.shape[0]/(pRbLo_mlcHi_Cells.shape[0] + pRbHi_mlcHi_Cells.shape[0] + pRbLo_mlcLo_Cells.shape[0] + pRbHi_mlcLo_Cells.shape[0]))
                    mlcLo.append(pRbHi_mlcLo_Cells.shape[0]/(pRbLo_mlcHi_Cells.shape[0] + pRbHi_mlcHi_Cells.shape[0] + pRbLo_mlcLo_Cells.shape[0] + pRbHi_mlcLo_Cells.shape[0]))






                    print(plate)
                    print(drug)
                    print(obs)
                    print(mlcHi)
                    print(mlcLo)



                    MLC2aHiCells = df.loc[((df['Mlc2a_Cyto'] >= mlc2a_gate) | (df['MLC2v_Cyto'] >= mlc2v_gate))]
                    MLC2aLoCells = df.loc[((df['Mlc2a_Cyto'] < mlc2a_gate) & (df['MLC2v_Cyto'] < mlc2v_gate))]
                    mlc_max = df[['Mlc2a_Cyto','MLC2v_Cyto']].max(axis=1)
                    df['mlc_max'] = mlc_max
                    df.loc[df['mlc_max'] == max(df['mlc_max'])] = np.NaN
                    plt=df.plot.scatter('Ki67_Nuc','mlc_max')
                    plt.axvline(ki76_gate)
                    plt.axhline(mlc2v_gate)
                    plt.set_title(drug.upper())
                    fig = plt.figure
                    fig.savefig('%s_%s_Ki67_%s.png' % (plate,obs,drug),dpi=600)

                    mlcHi = []
                    mlcLo = []
                    oct4aHi = []
                    oct4aLo = []

                    ki67Lo_mlcHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['Ki67_Nuc'] < ki76_gate]
                    ki67Hi_mlcHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['Ki67_Nuc'] >= ki76_gate]

                    ki67Lo_mlcLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['Ki67_Nuc'] < ki76_gate]
                    ki67Hi_mlcLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['Ki67_Nuc'] >= ki76_gate]

                    mlcHi.append(ki67Lo_mlcHi_Cells.shape[0]/(ki67Lo_mlcHi_Cells.shape[0] + ki67Hi_mlcHi_Cells.shape[0] + ki67Lo_mlcLo_Cells.shape[0] + ki67Hi_mlcLo_Cells.shape[0]))
                    mlcHi.append(ki67Hi_mlcHi_Cells.shape[0]/(ki67Lo_mlcHi_Cells.shape[0] + ki67Hi_mlcHi_Cells.shape[0] + ki67Lo_mlcLo_Cells.shape[0] + ki67Hi_mlcLo_Cells.shape[0]))

                    mlcLo.append(ki67Lo_mlcLo_Cells.shape[0]/(ki67Lo_mlcHi_Cells.shape[0] + ki67Hi_mlcHi_Cells.shape[0] + ki67Lo_mlcLo_Cells.shape[0] + ki67Hi_mlcLo_Cells.shape[0]))
                    mlcLo.append(ki67Hi_mlcLo_Cells.shape[0]/(ki67Lo_mlcHi_Cells.shape[0] + ki67Hi_mlcHi_Cells.shape[0] + ki67Lo_mlcLo_Cells.shape[0] + ki67Hi_mlcLo_Cells.shape[0]))


                    print(plate)
                    print(drug)
                    print(obs)
                    print(mlcHi)
                    print(mlcLo)








#                pRbLo_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] < pRb_gate]
#                pRbHi_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] >= pRb_gate]

#                pRbLo_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] < pRb_gate]
#                pRbHi_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] >= pRb_gate]

#                mlc2aHi.append(pRbLo_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
#                mlc2aHi.append(pRbHi_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

#                mlc2aLo.append(pRbLo_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
#                mlc2aLo.append(pRbHi_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

#                print(plate)
#                print(drug)
#                print(obs)
#                print(mlc2aHi)
#                print(mlc2aLo)






##Oct4a_Nuc
##Ki67_Nuc


#first = 1
#for plate in ['P6','P7','P8']:
#    for dose in ['lowdose','highdose']:
#        for ab in ['proteome1']:
#            df = pd.read_csv('%s/dmso/%s_dmso_ProcessedData%s_%s.csv' % (dose,dose,plate,ab))
#            if first == 1:  
#                pRb = df['Oct4a_Nuc']    #Try here 
#                mlc2v = df['Ki67_Nuc']
#                first = 0
#            else:
#                pRb = pRb.append(df['Oct4a_Nuc'])
#                mlc2v = mlc2v.append(df['Ki67_Nuc'])
#pRb = pRb.reset_index()
#mlc2v = mlc2v.reset_index()


#pRb.columns = ['index',0]
#pRb_threshold = gate_on_var(pRb)



#pRb_threshold = 12.5
##hplt=pRb[0].hist(bins=50)
##hplt.axvline(pRb_threshold)
##fig = hplt.figure
##fig.savefig(('Oct4a_Nuc.png'))
##matplotlib.pyplot.close()
##print(pRb_threshold)



#mlc2v.columns = ['index',0]
#mlc2v_threshold = gate_on_var(mlc2v)
#hplt2=mlc2v[0].hist(bins=25)
#hplt2.axvline(mlc2v_threshold)
#fig2 = hplt2.figure
#fig2.savefig(('Ki67_Nuc.png'))
#matplotlib.pyplot.close()



#    

