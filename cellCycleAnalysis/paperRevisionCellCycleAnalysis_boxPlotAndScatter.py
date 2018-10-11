Prot3
CNX43
pRb
MLC2v

Prot1 
a-ACTININ
MLC2v
p-Rb
Oct4a





import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


for plate in ['P6','P7','P8']:
    df = pd.read_csv('highdose/dmso/highdose_dmso_ProcessedData%s_proteome1.csv' % plate)
    plt=df.plot.scatter('Oct4a_Cyto','a-ACTININ_Cyto')
    fig = plt.figure
    fig.savefig(('%s_oct4a_aactinin.png' % plate))

    plt=df.plot.scatter('p-Rb_Cyto','a-ACTININ_Cyto')
    fig = plt.figure
    fig.savefig(('%s_p-Rb_aactinin.png' % plate))

    plt=df.plot.scatter('Oct4a_Cyto','MLC2v_Cyto')
    fig = plt.figure
    fig.savefig(('%s_oct4a_mlc2v.png' % plate))

    plt=df.plot.scatter('p-Rb_Cyto','MLC2v_Cyto')
    fig = plt.figure
    fig.savefig(('%s_p-Rb_mlc2v.png' % plate))


    plt=df.plot.scatter('p-Rb_Nuc','a-ACTININ_Cyto')
    fig = plt.figure
    fig.savefig(('%s_p-Rb_Nuc_aactinin.png' % plate))

    plt=df.plot.scatter('p-Rb_Nuc','MLC2v_Cyto')
    fig = plt.figure
    fig.savefig(('%s_p-Rb_Nuc_mlc2v.png' % plate))



for plate in ['P6','P7','P8']:
    df = pd.read_csv('highdose/dmso/highdose_dmso_ProcessedData%s_proteome3.csv' % plate)
    plt=df.plot.scatter('pRb_Cyto','CNX43_Cyto')
    fig = plt.figure
    fig.savefig(('%s_pRb_CNX43.png' % plate))

    plt=df.plot.scatter('pRb_Cyto','MLC2v_Cyto')
    fig = plt.figure
    fig.savefig(('%s_pRb_MLC2v.png' % plate))

    plt=df.plot.scatter('pRb_Nuc','CNX43_Cyto')
    fig = plt.figure
    fig.savefig(('%s_pRb_Nuc_CNX43.png' % plate))

    plt=df.plot.scatter('pRb_Nuc','MLC2v_Cyto')
    fig = plt.figure
    fig.savefig(('%s_pRb_Nuc_MLC2v.png' % plate))


#oct_cyt = df['Oct4a_Cyto']
#oct_nuc = df['Oct4a_Nuc']
#oct_combined = (oct_cyt+oct_nuc)/2



#act_cyt = df['a-ACTININ_Cyto']
#act_nuc = df['a-ACTININ_Nuc']
#act_combined = (act_cyt+act_nuc)/2

#split (cluster) on hi/lo cardiac marker 
#Time course of cell cycle marker in both conditions

#Prot3
#CNX43
#pRb
#MLC2v

#Prot1 
#a-ACTININ
#MLC2v
#p-Rb
#Oct4a

#for plate in ['P6','P7','P8']:
#    df = pd.read_csv('highdose/dmso/highdose_dmso_ProcessedData%s_proteome1.csv' % plate)


from sklearn.cluster import KMeans


df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP6_proteome1.csv')

#cluster on mlc2v
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_
print(centers)
#Identify hi and lo clusters
MLC2vHi_label = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
MLC2vLo_label = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])
print(MLC2vHi_label)


#Split into mlc2v_hi df
idxHi = clusters.index[clusters['cluster']==MLC2vHi_label].tolist()
MLC2vHiCells = df.loc[idxHi]
#And mlc2v_lo df
idxLo = clusters.index[clusters['cluster']==MLC2vLo_label].tolist()
MLC2vLoCells = df.loc[idxLo]

print(MLC2vHiCells.shape[0])
print(MLC2vLoCells.shape[0])



#re-gate MLC2v hi and lo on pRb
#First MLC2vHi 
model = KMeans(n_clusters=2,init='k-means++')
model.fit(MLC2vHiCells['p-Rb_Nuc'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
pRbHi_MLC2vHi_Cells = df.loc[idxHi].shape[0]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
pRbLo_MLC2vHi_Cells = df.loc[idxLo].shape[0]

#Then MLC2vLo 

model = KMeans(n_clusters=2,init='k-means++')
model.fit(MLC2vLoCells['p-Rb_Nuc'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
pRbHi_MLC2vLo_Cells = df.loc[idxHi].shape[0]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
pRbLo_MLC2vLo_Cells = df.loc[idxLo].shape[0]
#print(pRbHiCells.shape[0])

pRb_hi_mlc2v_hi_percent = pRbHi_MLC2vHi_Cells/(pRbHi_MLC2vHi_Cells+pRbLo_MLC2vHi_Cells)
pRb_hi_mlc2v_lo_percent = pRbHi_MLC2vLo_Cells/(pRbHi_MLC2vLo_Cells+pRbLo_MLC2vLo_Cells)

print(pRb_hi_mlc2v_hi_percent)
print(pRb_hi_mlc2v_lo_percent)





pRb_hi_oct4a_hi_p6 = Oct4aHiCells['p-Rb_Nuc']
pRb_lo_p6 = Oct4aLoCells['p-Rb_Nuc']


pRb_lo_p6 = Oct4aLoCells['p-Rb_Nuc']

df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP7_proteome1.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['a-ACTININ_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p7 = Oct4aHiCells['p-Rb_Nuc']
pRb_lo_p7 = Oct4aLoCells['p-Rb_Nuc']



#dfNew = pd.concat([pRb_hi_p6,pRb_lo_p6], ignore_index=True, axis=1)
#dfNew.columns = ['a-ACTININ_HiCells_P6','a-ACTININ_LoCells_P6']

df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP8_proteome1.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['a-ACTININ_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p8 = Oct4aHiCells['p-Rb_Nuc']
pRb_lo_p8 = Oct4aLoCells['p-Rb_Nuc']


dfNew = pd.concat([pRb_hi_p6,pRb_lo_p6,pRb_hi_p7,pRb_lo_p7,pRb_hi_p8,pRb_lo_p8], ignore_index=True, axis=1)
dfNew.columns = ['ACTN1_hi_P6','ACTN1_lo_P6','ACTN1_hi_P7','ACTN1_lo_P7','ACTN1_hi_P8','ACTN1_lo_P8']


plt=dfNew.plot.box()
plt.set_title = 'pRb'
fig = plt.figure
fig.savefig('pRb_actn1_boxplots_soraf.png')

#a-ACTININ
#MLC2v
#p-Rb
#Oct4a

###################################################



from sklearn.cluster import KMeans
df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP1_signaling1.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p6 = Oct4aHiCells['p-Rb_Nuc']
pRb_lo_p6 = Oct4aLoCells['p-Rb_Nuc']




df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP2_signaling1.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p7 = Oct4aHiCells['p-Rb_Nuc']
pRb_lo_p7 = Oct4aLoCells['p-Rb_Nuc']



#dfNew = pd.concat([pRb_hi_p6,pRb_lo_p6], ignore_index=True, axis=1)
#dfNew.columns = ['a-ACTININ_HiCells_P6','a-ACTININ_LoCells_P6']

df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP4_signaling1.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p8 = Oct4aHiCells['p-Rb_Nuc']
pRb_lo_p8 = Oct4aLoCells['p-Rb_Nuc']


dfNew = pd.concat([pRb_hi_p6,pRb_lo_p6,pRb_hi_p7,pRb_lo_p7,pRb_hi_p8,pRb_lo_p8], ignore_index=True, axis=1)
dfNew.columns = ['MLC2v_hi_P1','MLC2v_lo_P1','MLC2v_hi_P2','MLC2v_lo_P2','MLC2v_hi_P4','MLC2v_lo_P4']


plt=dfNew.plot.box()
plt.set_title = 'pRb'
fig = plt.figure
fig.savefig('pRb_MLC2v_boxplots_soraf.png')



###############################################



from sklearn.cluster import KMeans
df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP1_signaling2.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p6 = Oct4aHiCells['Oct4a_Nuc']
pRb_lo_p6 = Oct4aLoCells['Oct4a_Nuc']




df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP2_signaling2.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p7 = Oct4aHiCells['Oct4a_Nuc']
pRb_lo_p7 = Oct4aLoCells['Oct4a_Nuc']



#dfNew = pd.concat([pRb_hi_p6,pRb_lo_p6], ignore_index=True, axis=1)
#dfNew.columns = ['a-ACTININ_HiCells_P6','a-ACTININ_LoCells_P6']

df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP4_signaling2.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p8 = Oct4aHiCells['Oct4a_Nuc']
pRb_lo_p8 = Oct4aLoCells['Oct4a_Nuc']


dfNew = pd.concat([pRb_hi_p6,pRb_lo_p6,pRb_hi_p7,pRb_lo_p7,pRb_hi_p8,pRb_lo_p8], ignore_index=True, axis=1)
dfNew.columns = ['MLC2v_hi_P1','MLC2v_lo_P1','MLC2v_hi_P2','MLC2v_lo_P2','MLC2v_hi_P4','MLC2v_lo_P4']


plt=dfNew.plot.box()
plt.set_title = 'Oct4a_Nuc'
fig = plt.figure
fig.savefig('oct4a_MLC2v_boxplots_soraf.png')





##############################

from sklearn.cluster import KMeans
df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP1_signaling2.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p6 = Oct4aHiCells['pRb_Nuc']
pRb_lo_p6 = Oct4aLoCells['pRb_Nuc']




df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP2_signaling2.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p7 = Oct4aHiCells['pRb_Nuc']
pRb_lo_p7 = Oct4aLoCells['pRb_Nuc']



#dfNew = pd.concat([pRb_hi_p6,pRb_lo_p6], ignore_index=True, axis=1)
#dfNew.columns = ['a-ACTININ_HiCells_P6','a-ACTININ_LoCells_P6']

df = pd.read_csv('highdose/soraf/highdose_soraf_ProcessedDataP4_signaling2.csv')
#Maybe do this with nuc pRb col in whole df, then plot box plot of results
model = KMeans(n_clusters=2,init='k-means++')
model.fit(df['MLC2v_Cyto'].values.reshape(-1,1))
labels = model.labels_  
clusters = pd.DataFrame(data=labels,columns=['cluster'])
centers = model.cluster_centers_

pRbHi = np.argmax([np.mean(centers[0,]),np.mean(centers[1,])])
pRbLo = np.argmin([np.mean(centers[0,]),np.mean(centers[1,])])

idxHi = clusters.index[clusters['cluster']==pRbHi].tolist()
Oct4aHiCells = df.loc[idxHi]

idxLo = clusters.index[clusters['cluster']==pRbLo].tolist()
Oct4aLoCells = df.loc[idxLo]


pRb_hi_p8 = Oct4aHiCells['pRb_Nuc']
pRb_lo_p8 = Oct4aLoCells['pRb_Nuc']


dfNew = pd.concat([pRb_hi_p6,pRb_lo_p6,pRb_hi_p7,pRb_lo_p7,pRb_hi_p8,pRb_lo_p8], ignore_index=True, axis=1)
dfNew.columns = ['MLC2v_hi_P1','MLC2v_lo_P1','MLC2v_hi_P2','MLC2v_lo_P2','MLC2v_hi_P4','MLC2v_lo_P4']


plt=dfNew.plot.box()
plt.set_title = 'pRb'
fig = plt.figure
fig.savefig('pRb_MLC2v_boxplots_soraf_sig2.png')




