from sklearn.cluster import KMeans
from statistics import mode
from statistics import median
from statistics import mean
from statistics import StatisticsError
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#one line - plate 1,2,4
#Sig1  MLC2v_Cyto,p_Rb_Nuc
#Sig2  MLC2v_Cyto, pRB_Nuc
#Dose - low, medium, high 

#Combine all data, cluster 2 var, get one threshold for gating.
#Use DMSO for gating
#Then plot percents 2 panels x 3 doses = 6 lines (6 rows used on plates 1-4)
#First DMSO, then soraf 

#    df = pd.read_csv('highdose/dmso/highdose_dmso_ProcessedData%s_proteome1.csv' % plate)
#    plt=df.plot.scatter('Oct4a_Cyto','a-ACTININ_Cyto')
#    fig = plt.figure
#    fig.savefig(('%s_oct4a_aactinin.png' % plate))
#Mlc2a, proteome1 and 2

#def gate_on_var(df,var):
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


        percentHi = len(idxHi)/(len(idxLo)+len(idxHi))
        percentHi_all.append(percentHi)

        try:
            percentHi_final = mode(percentHi_all)
        except StatisticsError:
            percentHi_final = median(percentHi_all)

    threshold = mean([min(hiCells[0]),max(loCells[0])])

    return threshold




first = 1
for plate in ['P1','P2','P4']:
    for dose in ['lowdose','meddose','highdose']:
        for ab in ['signaling1','signaling2']:
            df = pd.read_csv('%s/dmso/%s_dmso_ProcessedData%s_%s.csv' % (dose,dose,plate,ab))
            if first == 1:  
                pRb = df['p-Rb_Nuc']    #Try here 
                mlc2v = df['MLC2v_Cyto']
                first = 0
            else:
                try:
                    pRb = pRb.append(df['p-Rb_Nuc'])

                except:
                    pRb = pRb.append(df['pRb_Nuc'])
                mlc2v = mlc2v.append(df['MLC2v_Cyto'])


hplt=pRb.hist(bins=25)
fig = hplt.figure
fig.savefig(('pRbHist.png'))

hplt2=mlc2v.hist(bins=25)
fig2 = hplt2.figure
fig2.savefig(('mlc2vHist.png'))

pRb = pRb.reset_index()
mlc2v = mlc2v.reset_index()
pRb_threshold = gate_on_var(pRb)
mlc2v.columns = ['index',0]
mlc2v_threshold = gate_on_var(mlc2v)

#In [3]: pRb_threshold
#Out[3]: 10.408327849108183

#In [4]: mlc2v_threshold
#Out[4]: 11.554632675741644







###########################################3
#Proteome





for dose in ['lowdose','highdose']:
    for ab in ['proteome1','proteome3']:
        three_time_points_cardiac = []
        three_time_points_noncardiac = []
        three_time_points_control = []
        for plate in ['P6','P7','P8']:
            df = pd.read_csv('%s/soraf/%s_soraf_ProcessedData%s_%s.csv' % (dose,dose,plate,ab))

            if df.shape[0] > 0: #lowdose_soraf_ProcessedDataP2_signaling2.csv is empty

                MLC2vHiCells = df.loc[df['MLC2v_Cyto'] >= mlc2v_threshold]
                MLC2vLoCells = df.loc[df['MLC2v_Cyto'] < mlc2v_threshold]
                #percentHi_cardiac = pRbHi_MLC2vHi_Cells.shape[0]/(pRbHi_MLC2vHi_Cells.shape[0] + pRbHi_MLC2vLo_Cells.shape[0])

                try:    
                    pRbHi_MLC2vHi_Cells = MLC2vHiCells.loc[MLC2vHiCells['p-Rb_Nuc'] >= pRb_threshold]
                    pRbLo_MLC2vHi_Cells = MLC2vHiCells.loc[MLC2vHiCells['p-Rb_Nuc'] < pRb_threshold]
                    percentHi_cardiac = pRbHi_MLC2vHi_Cells.shape[0]/(pRbHi_MLC2vHi_Cells.shape[0] + pRbLo_MLC2vHi_Cells.shape[0])

                    pRbHi_MLC2vLo_Cells = MLC2vLoCells.loc[MLC2vLoCells['p-Rb_Nuc'] >= pRb_threshold]
                    pRbLo_MLC2vLo_Cells = MLC2vLoCells.loc[MLC2vLoCells['p-Rb_Nuc'] < pRb_threshold]
                    percentHi_noncardiac = pRbHi_MLC2vLo_Cells.shape[0]/(pRbHi_MLC2vLo_Cells.shape[0] + pRbLo_MLC2vLo_Cells.shape[0])

                    pRbHi = df.loc[df['p-Rb_Nuc'] >= pRb_threshold]
                    pRbLo = df.loc[df['p-Rb_Nuc'] < pRb_threshold]
                    percentHi_all = pRbHi.shape[0]/(pRbHi.shape[0] + pRbLo.shape[0])

                except KeyError:
                    pRbHi_MLC2vHi_Cells = MLC2vHiCells.loc[MLC2vHiCells['pRb_Nuc'] >= pRb_threshold]
                    pRbLo_MLC2vHi_Cells = MLC2vHiCells.loc[MLC2vHiCells['pRb_Nuc'] < pRb_threshold]
                    percentHi_cardiac = pRbHi_MLC2vHi_Cells.shape[0]/(pRbHi_MLC2vHi_Cells.shape[0] + pRbHi_MLC2vLo_Cells.shape[0])

                    pRbHi_MLC2vLo_Cells = MLC2vLoCells.loc[MLC2vLoCells['pRb_Nuc'] >= pRb_threshold]
                    pRbLo_MLC2vLo_Cells = MLC2vLoCells.loc[MLC2vLoCells['pRb_Nuc'] < pRb_threshold]
                    percentHi_noncardiac = pRbHi_MLC2vLo_Cells.shape[0]/(pRbHi_MLC2vLo_Cells.shape[0] + pRbLo_MLC2vLo_Cells.shape[0])

                    pRbHi = df.loc[df['pRb_Nuc'] >= pRb_threshold]
                    pRbLo = df.loc[df['pRb_Nuc'] < pRb_threshold]
                    percentHi_all = pRbHi.shape[0]/(pRbHi.shape[0] + pRbLo.shape[0])

                three_time_points_cardiac.append(percentHi_cardiac)
                three_time_points_noncardiac.append(percentHi_noncardiac)
                three_time_points_control.append(percentHi_all)

            else:
                three_time_points_cardiac.append(0)
                three_time_points_noncardiac.append(0)
                three_time_points_control.append(0)

        fig = plt.figure()
        x = [1,2,24]
        plt.plot(x,three_time_points_cardiac,'g',label='MLC2v+')
        plt.plot(x,three_time_points_noncardiac,'r',label='MLC2v-')
        plt.plot(x,three_time_points_control,'b')

        plt.xlabel('Time hrs')
        plt.ylabel('% pRb Positive Cells')
        plt.savefig('%s_%s.png' % (dose,ab))



for dose in ['lowdose','highdose']:
    for ab in ['proteome1','proteome3']:
        three_time_points_cardiac = []
        three_time_points_noncardiac = []
        three_time_points_control = []
        for plate in ['P6','P7','P8']:
            df = pd.read_csv('%s/dmso/%s_dmso_ProcessedData%s_%s.csv' % (dose,dose,plate,ab))

            if df.shape[0] > 0: #lowdose_soraf_ProcessedDataP2_signaling2.csv is empty

                MLC2vHiCells = df.loc[df['MLC2v_Cyto'] >= mlc2v_threshold]
                MLC2vLoCells = df.loc[df['MLC2v_Cyto'] < mlc2v_threshold]
                #percentHi_cardiac = pRbHi_MLC2vHi_Cells.shape[0]/(pRbHi_MLC2vHi_Cells.shape[0] + pRbHi_MLC2vLo_Cells.shape[0])

                try:    
                    pRbHi_MLC2vHi_Cells = MLC2vHiCells.loc[MLC2vHiCells['p-Rb_Nuc'] >= pRb_threshold]
                    pRbLo_MLC2vHi_Cells = MLC2vHiCells.loc[MLC2vHiCells['p-Rb_Nuc'] < pRb_threshold]
                    percentHi_cardiac = pRbHi_MLC2vHi_Cells.shape[0]/(pRbHi_MLC2vHi_Cells.shape[0] + pRbLo_MLC2vHi_Cells.shape[0])

                    pRbHi_MLC2vLo_Cells = MLC2vLoCells.loc[MLC2vLoCells['p-Rb_Nuc'] >= pRb_threshold]
                    pRbLo_MLC2vLo_Cells = MLC2vLoCells.loc[MLC2vLoCells['p-Rb_Nuc'] < pRb_threshold]
                    percentHi_noncardiac = pRbHi_MLC2vLo_Cells.shape[0]/(pRbHi_MLC2vLo_Cells.shape[0] + pRbLo_MLC2vLo_Cells.shape[0])

                    pRbHi = df.loc[df['p-Rb_Nuc'] >= pRb_threshold]
                    pRbLo = df.loc[df['p-Rb_Nuc'] < pRb_threshold]
                    percentHi_all = pRbHi.shape[0]/(pRbHi.shape[0] + pRbLo.shape[0])

                except KeyError:
                    pRbHi_MLC2vHi_Cells = MLC2vHiCells.loc[MLC2vHiCells['pRb_Nuc'] >= pRb_threshold]
                    pRbLo_MLC2vHi_Cells = MLC2vHiCells.loc[MLC2vHiCells['pRb_Nuc'] < pRb_threshold]
                    percentHi_cardiac = pRbHi_MLC2vHi_Cells.shape[0]/(pRbHi_MLC2vHi_Cells.shape[0] + pRbHi_MLC2vLo_Cells.shape[0])

                    pRbHi_MLC2vLo_Cells = MLC2vLoCells.loc[MLC2vLoCells['pRb_Nuc'] >= pRb_threshold]
                    pRbLo_MLC2vLo_Cells = MLC2vLoCells.loc[MLC2vLoCells['pRb_Nuc'] < pRb_threshold]
                    percentHi_noncardiac = pRbHi_MLC2vLo_Cells.shape[0]/(pRbHi_MLC2vLo_Cells.shape[0] + pRbLo_MLC2vLo_Cells.shape[0])

                    pRbHi = df.loc[df['pRb_Nuc'] >= pRb_threshold]
                    pRbLo = df.loc[df['pRb_Nuc'] < pRb_threshold]
                    percentHi_all = pRbHi.shape[0]/(pRbHi.shape[0] + pRbLo.shape[0])

                three_time_points_cardiac.append(percentHi_cardiac)
                three_time_points_noncardiac.append(percentHi_noncardiac)
                three_time_points_control.append(percentHi_all)

            else:
                three_time_points_cardiac.append(0)
                three_time_points_noncardiac.append(0)
                three_time_points_control.append(0)

        fig = plt.figure()
        x = [1,2,24]
        plt.plot(x,three_time_points_cardiac,'g',label='MLC2v+')
        plt.plot(x,three_time_points_noncardiac,'r',label='MLC2v-')
        plt.plot(x,three_time_points_control,'b')

        plt.xlabel('Time hrs')
        plt.ylabel('% pRb Positive Cells')
        plt.savefig('%s_%s.png' % (dose,ab))





####################33
#Scatter plots

mlc2a_gate = 11.2
pRb_gate = 10.4
mlc2v_gate = 10.9
#plt.axvline(x_position)


for dose in ['lowdose','highdose']:
    for ab in ['proteome1']:
        three_time_points_cardiac = []
        three_time_points_noncardiac = []
        three_time_points_control = []
        for plate in ['P6','P7','P8']:
            mlc2aHi = []
            mlc2aLo = []

            df = pd.read_csv('%s/dmso/%s_dmso_ProcessedData%s_%s.csv' % (dose,dose,plate,ab))

            #df = pd.read_csv('highdose/dmso/highdose_dmso_ProcessedData%s_proteome1.csv' % plate)
            plt=df.plot.scatter('p-Rb_Nuc','Mlc2a_Cyto')
            plt.axvline(pRb_gate)
            plt.axhline(mlc2a_gate)
            fig = plt.figure
            fig.savefig('%s_%s_p-Rb_mlc2a.png' % (dose,plate))




            MLC2aHiCells = df.loc[df['Mlc2a_Cyto'] >= mlc2a_gate]
            MLC2aLoCells = df.loc[df['Mlc2a_Cyto'] < mlc2a_gate]


            pRbLo_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] >= pRb_gate]

            pRbLo_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] >= pRb_gate]

            mlc2aHi.append(pRbLo_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aHi.append(pRbHi_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            mlc2aLo.append(pRbLo_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aLo.append(pRbHi_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            print(plate)
            print(dose)
            print(mlc2aHi)
            print(mlc2aLo)


mlc2a_gate = 11.2
pRb_gate = 10.4
mlc2v_gate = 10.9
#plt.axvline(x_position)


for dose in ['lowdose','highdose']:
    for ab in ['proteome1']:
        three_time_points_cardiac = []
        three_time_points_noncardiac = []
        three_time_points_control = []
        for plate in ['P6','P7','P8']:
            mlc2aHi = []
            mlc2aLo = []

            df = pd.read_csv('%s/soraf/%s_soraf_ProcessedData%s_%s.csv' % (dose,dose,plate,ab))
            plt=df.plot.scatter('p-Rb_Nuc','Mlc2a_Cyto')
            plt.axvline(pRb_gate)
            plt.axhline(mlc2a_gate)
            fig = plt.figure
            fig.savefig('%s_%s_p-Rb_mlc2a.png' % (dose,plate))




            MLC2aHiCells = df.loc[df['Mlc2a_Cyto'] >= mlc2a_gate]
            MLC2aLoCells = df.loc[df['Mlc2a_Cyto'] < mlc2a_gate]


            pRbLo_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] >= pRb_gate]

            pRbLo_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] >= pRb_gate]

            mlc2aHi.append(pRbLo_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aHi.append(pRbHi_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            mlc2aLo.append(pRbLo_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aLo.append(pRbHi_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            print(plate)
            print(dose)
            print(mlc2aHi)
            print(mlc2aLo)



#Mlc2a, proteome1 and 2
#Oct4a prot1
#p-Rb prot1
##############################3
#Oct4a

#Get oct4a threshold
first = 1
for plate in ['P6','P7','P8']:
    for dose in ['lowdose','highdose']:
        for ab in ['proteome1']:
            df = pd.read_csv('%s/dmso/%s_dmso_ProcessedData%s_%s.csv' % (dose,dose,plate,ab))
            if first == 1:  
                oct4a = df['Oct4a_Nuc']
                first = 0
            else:
                oct4a = oct4a.append(df['Oct4a_Nuc'])


hplt=oct4a.hist(bins=25)
fig = hplt.figure
fig.savefig(('oct4aHist.png'))


oct4a = oct4a.reset_index()
oct4a.columns = ['index',0]
oct4a_gate = gate_on_var(oct4a)


oct4a_gate=11


for dose in ['lowdose','highdose']:
    for ab in ['proteome1']:
        three_time_points_cardiac = []
        three_time_points_noncardiac = []
        three_time_points_control = []
        for plate in ['P6','P7','P8']:
            mlc2aHi = []
            mlc2aLo = []

            df = pd.read_csv('%s/dmso/%s_dmso_ProcessedData%s_%s.csv' % (dose,dose,plate,ab))

            #df = pd.read_csv('highdose/dmso/highdose_dmso_ProcessedData%s_proteome1.csv' % plate)
            plt=df.plot.scatter('p-Rb_Nuc','Oct4a_Nuc')
            plt.axvline(pRb_gate)
            plt.axhline(oct4a_gate)
            fig = plt.figure
            fig.savefig('%s_%s_p-Rb_oct4a.png' % (dose,plate))




            MLC2aHiCells = df.loc[df['Oct4a_Nuc'] >= mlc2a_gate]
            MLC2aLoCells = df.loc[df['Oct4a_Nuc'] < mlc2a_gate]


            pRbLo_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] >= pRb_gate]

            pRbLo_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] >= pRb_gate]

            mlc2aHi.append(pRbLo_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aHi.append(pRbHi_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            mlc2aLo.append(pRbLo_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aLo.append(pRbHi_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            print(plate)
            print(dose)
            print(mlc2aHi)
            print(mlc2aLo)



for dose in ['lowdose','highdose']:
    for ab in ['proteome1']:
        three_time_points_cardiac = []
        three_time_points_noncardiac = []
        three_time_points_control = []
        for plate in ['P6','P7','P8']:
            mlc2aHi = []
            mlc2aLo = []

            df = pd.read_csv('%s/soraf/%s_soraf_ProcessedData%s_%s.csv' % (dose,dose,plate,ab))
            plt=df.plot.scatter('p-Rb_Nuc','Oct4a_Nuc')
            plt.axvline(pRb_gate)
            plt.axhline(oct4a_gate)
            fig = plt.figure
            fig.savefig('%s_%s_p-Rb_oct4a.png' % (dose,plate))




            MLC2aHiCells = df.loc[df['Oct4a_Nuc'] >= mlc2a_gate]
            MLC2aLoCells = df.loc[df['Oct4a_Nuc'] < mlc2a_gate]


            pRbLo_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aHi_Cells = MLC2aHiCells.loc[MLC2aHiCells['p-Rb_Nuc'] >= pRb_gate]

            pRbLo_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] < pRb_gate]
            pRbHi_MLC2aLo_Cells = MLC2aLoCells.loc[MLC2aLoCells['p-Rb_Nuc'] >= pRb_gate]

            mlc2aHi.append(pRbLo_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aHi.append(pRbHi_MLC2aHi_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            mlc2aLo.append(pRbLo_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))
            mlc2aLo.append(pRbHi_MLC2aLo_Cells.shape[0]/(pRbHi_MLC2aHi_Cells.shape[0] + pRbHi_MLC2aLo_Cells.shape[0] + pRbLo_MLC2aHi_Cells.shape[0] + pRbLo_MLC2aLo_Cells.shape[0]))

            print(plate)
            print(dose)
            print(mlc2aHi)
            print(mlc2aLo)

