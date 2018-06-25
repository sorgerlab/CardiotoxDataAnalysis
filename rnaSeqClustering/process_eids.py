import glob

#Matlab only processes and ouputs the numbers at tail end of ensemble id 
#Append 'ensg' tag and leading 0's to allow further processing 
for filename in sorted(glob.iglob('ensembleids/*.csv')):
    #print('%s' % filename)
    newlist = []
    f = open('%s' % filename,'r')

    for line in f:
        line = line.rstrip()
        line = line.split('.')[0]
        toadd = 11-len(line)
        newline = 'ENSG' +'0'*toadd + line
        newlist.append(newline)
    f.close()
    index = filename.split('_')[2].split('.')[0]
    print(index)
    with open('ensembleids_processed/genes_cluster_' '%s' % index + '.csv' ,'w') as newfile:
        newfile.write('\n'.join(newlist))



