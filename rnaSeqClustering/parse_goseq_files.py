import glob

#Loop through all files generated in goseq analysis. 
#Searching for top go terms and associated genes
for filename in sorted(glob.iglob('*goseq.txt')):
    gene_number = []
    genes = []
    goid = []
    gonames = []
    with open(filename,'r') as f:
        goseqdata = f.readlines()
        tophit = goseqdata[1].split('\t')
        secondhit = goseqdata[2].split('\t')
        thirdhit = goseqdata[3].split('\t')
        goid.append(tophit[0])
        gonames.append(tophit[1])
        genes = tophit[6].strip().split(', ') #if reading from file
        #Looking for 20 genes to plot, if less than 20 genes in top go term, proceed to next
        if len(genes) < 20:
            genes = genes + secondhit[6].strip().split(', ')
            goid.append(secondhit[0] + '*')
            gonames.append(secondhit[1] + '*')
        if len(genes) < 20:
            genes = genes + thirdhit[6].strip().split(', ')
            goid.append(thirdhit[0] + '*')   
            gonames.append(thirdhit[1] + '*')     
        genes = genes[0:20]
        gene_number.append(len(genes))
        print(filename)

    #Save ensemble ids for genes, associated go id number and go term name
    with open('ensemble_ids.txt','a') as f2:
        f2.write('\n'.join(genes))
        f2.write('\n')
    with open('go_ids.txt','a') as f2:
        f2.write('\n'.join(goid))
        f2.write('\n')
    with open('go_terms.txt','a') as f2:
        f2.write('\n'.join(gonames))
        f2.write('\n')


##

#Take identified ensemble ids for each cluster 
#Find the index for this gene in matrix of all differentially expressed genes 
#Use indices for plotting in matlab 
with open('ensembleids_de_all.csv','r') as f:
    all_de_genes = f.read().splitlines()
clusterGeneIndex = []
with open('ensemble_ids_clusters.txt','r') as fgenes:
    genes = fgenes.read().splitlines()
for gene in genes:
    if gene in all_de_genes:
        clusterGeneIndex.append(all_de_genes.index(gene)+1)
with open('clusterGeneIndex.csv','w') as f2:
    for ind in clusterGeneIndex:
        f2.write("%i\n" % ind)


