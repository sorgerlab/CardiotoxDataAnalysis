import pandas as pd

df1 = pd.read_csv('rae_Cardio_Day1_5treatments_ProteinQuant_noscale_upid_ibaq_nologtrans_scaled.csv')
df2 = pd.read_csv('rae_Cardio_Day3_5treatments_ProteinQuant_noscale_upid_ibaq_nologtrans_scaled.csv')
#genes1 = list(df1['Gene Symbol']) 
#genes2 = list(df2['Gene Symbol']) 

#union = [el for el in genes1 if el in genes2]                      
#union = [el for el in union if el in genes1]

#dfunion1 = df1[df1['Gene Symbol'].isin(union)]
#dfunion2 = df2[df2['Gene Symbol'].isin(union)]

#df1_unqique = dfunion1.drop_duplicates(['Gene Symbol'])
#df2_unqique = dfunion2.drop_duplicates(['Gene Symbol'])


#df1_unqique.to_csv('day1_ibaq_scaled_matching.csv')
#df2_unqique.to_csv('day3_ibaq_scaled_matching.csv')


genes1 = list(df1['Uniprot_Id']) 
genes2 = list(df2['Uniprot_Id']) 

union = [el for el in genes1 if el in genes2]                      
union = [el for el in union if el in genes1]

dfunion1 = df1[df1['Uniprot_Id'].isin(union)]
dfunion2 = df2[df2['Uniprot_Id'].isin(union)]

df1_unqique = dfunion1.drop_duplicates(['Uniprot_Id'])
df2_unqique = dfunion2.drop_duplicates(['Uniprot_Id'])


df1_unqique.to_csv('day1_ibaq_scaled_matching.csv')
df2_unqique.to_csv('day3_ibaq_scaled_matching.csv')

