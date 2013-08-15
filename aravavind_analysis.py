# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import synapseclient
import IPython.display as display
import numpy as np
import pandas as pd
import scipy.stats as stats
from MicroArray import QaD_SVD, scale

DATA_ID = 'syn1968267'
METADATA_ID = 'syn2024470'
syn = synapseclient.login()

# <markdowncell>

# #Download and annotate data
# 1. Downloads expression data from Synapse
# 2. Gets annotations from metadata in Synapse
# 3. Extracts the data that we want to analyze
# 4. log transform data

# <codecell>

dataEnt=syn.get(DATA_ID)
#Load the data
df=pd.read_csv(dataEnt.path,sep='\t', index_col=0)
genes, geneNames, geneLoci = np.asarray(df.index, dtype=np.str), np.asarray(df['symbol'], dtype=np.str), np.asarray(df['locus'], dtype=np.str)
df = df.drop(['symbol', 'locus'], axis=1)
labels = np.asarray(df.keys())
data = df.as_matrix()

# <markdowncell>

# Query synapse for metadata of the expression data and store in a pandas dataframe

# <codecell>

#Specifically select those samples that are Stem Cells induced by "OSKM-NLT" genes
metadata = syn.query('select * from entity where project=="Arv" and parentId=="%s"' %METADATA_ID)['results']
#Add the hESC SC samples
#metadata.extend(syn.query('select * from entity where diffnameshort=="SC" and celltypeoforigin=="hESC" and parentId=="%s"' %METADATA_ID)['results'])
metadata = pd.DataFrame(metadata)
#Remove the unecessary lists and 'entity' in names  (this should be fixed on Synapse!)
for key in metadata.keys():
    if key in ['entity.benefactorId', 'entity.concreteType', 'entity.createdByPrincipalId', 'entity.createdOn', 'entity.createdByPrincipalId', 'entity.id', 'entity.modifiedOn', 'entity.modifiedByPrincipalId', 'entity.noteType', 'entity.versionLabel', 'entity.versionComment', 'entity.versionNumber', 'entity.parentId', 'entity.description']:
        del metadata[key]
        continue
    newkey=key.replace('entity.', '')
    metadata[newkey] = pd.Series([item[0] for item in metadata[key] if type(item) is list])
    del metadata[key]

#Display a table of metadata
pd.set_printoptions(max_colwidth=-1)
display.HTML(metadata[['bamName', 'celltypeoforigin', 'celltissuesourcepreinduction', 'disease', 'diffnameshort', 'run', 'lane', 'index', 'donorsex.cell.lines', 'numberofreads', 'uniquelyalignedreads', 'ratio']].to_html())

# <markdowncell>

# Filter out the Samples corresponding to found metadata

# <codecell>

idx = np.where([label in list(metadata['bamName']) for label in labels])[0]
labels = labels[idx]
data = data[:, idx]
print data.shape

# <markdowncell>

# Fillter out those genes where more than 90% of the samples have 0 expression and variance <=.5
# Then normalize by subtracting the mean of each row.
# Visualize

# <codecell>

#Filter by:        number of 0 items for each gene             variance being larger than...
idx=np.logical_and(np.sum(data==0, 1)/float(data.shape[1])<.9, np.var(np.log2(data+1),1)>.2)
data=data[idx]
genes = genes[idx]
geneNames = geneNames[idx]
geneLoci = geneLoci[idx]

#Log transform and mean Center
dataNorm=np.log2(data+1)
dataNorm=scale(dataNorm, center=True, scale=True)
print dataNorm.shape

#Plot distribution and standard deviation
mu, sigma=stats.norm.fit(dataNorm.flatten())
pylab.figure(figsize=(10,3));
pylab.subplot(1,2,1);
n, bins, patches = pylab.hist(dataNorm.flatten(), 100, normed=True);
pylab.plot(bins, stats.norm.pdf(bins, mu, sigma), 'r', linewidth=2)
pylab.subplot(1,2,2);
pylab.hist(np.std(dataNorm,1),100);

# <codecell>

metaColor = metadata['celltissuesourcepreinduction']
u, s, vt = QaD_SVD(dataNorm, colorLabels = metaColor, labels=metadata['bamName'])

# <codecell>


# <markdowncell>

# Output the top associated genes with each 

# <codecell>

topGenes=pd.DataFrame(range(100), columns=['none'])
topGeneLabels=[]
for i in range(3):
    topGeneIdx=argsort(u[:,i])
    topGenes['AntiCorr Eigengene %i' %(i+1)] = geneNames[topGeneIdx][:100]
    topGenes['Correlated Eigengene %i' %(i+1)] = geneNames[topGeneIdx][-100:]
del topGenes['none']
display.HTML(topGenes.to_html())


