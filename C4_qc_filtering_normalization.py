# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import synapseclient
import numpy as np
import pandas as pd
import scipy.stats as stats
import hierarchical_clustering as clust
import MicroArray
from toppGenePost import ToppGeneEnrichement
from IPython.display import HTML

EXPR_ID = 'syn2247799' #private'syn1968267'
EXPR_META_ID = 'syn2278178'
METH_ID = 'syn2233188'
METH_META_ID = 'syn2677043'
MIRSEQ_ID = 'syn2247832' #Private 'syn2233189'
MIRSEQ_META_ID = 'syn2278179'
BIOMART_ANNOT_ID = 'syn2706233'
pd.set_option('max_columns', 200)
syn = synapseclient.login()

# <markdowncell>

# #Download and annotate data
# 1. Downloads expression/miRNA-Seq/methylation data from Synapse
# 2. Gets annotations from metadata in Synapse
# 3. Make sure the metadata matches up with samples

# <codecell>

mirseq = pd.read_csv(syn.get(MIRSEQ_ID).path, sep='\t', index_col=0)
mirseq_meta = pd.read_csv(syn.get(MIRSEQ_META_ID).path, sep='\t',index_col=0)
mirseq_meta = mirseq_meta.ix[mirseq.columns,:]
mirseq_meta.drop(['Related C4 ID comment'],1, inplace=True)
print 'Mirseq Data Size:', mirseq.shape, mirseq_meta.shape
rnaseq = pd.read_csv(syn.get(EXPR_ID).path, sep='\t', index_col=[0,1,2])
rnaseq_meta = pd.read_csv(syn.get(EXPR_META_ID).path, sep='\t',index_col=0)
rnaseq_meta = rnaseq_meta.ix[rnaseq.columns,:]
print 'RNA-Seq Data Size:', rnaseq.shape, rnaseq_meta.shape
methyl = pd.read_csv(syn.get(METH_ID).path, sep='\t', index_col=0)
methyl_meta = pd.read_csv(syn.get(METH_META_ID).path, sep='\t',index_col=0)
methyl_meta = methyl_meta.ix[methyl.columns,:]
print 'Methylation Data Size:', methyl.shape, methyl_meta.shape

# <markdowncell>

# #Data exploration, Quality assesement and Normalization

# <markdowncell>

# ###micro RNA-Seq

# <codecell>

mirseq[mirseq.isnull()]=0
#microRNA probes where >10% of samples have more than 2 reads
idx = (sum(mirseq<2, axis=1)/mirseq.shape[1])<=.90
mirseq_norm=mirseq[idx]
print 'filterted down to', sum(idx), 'probes based on low coverage probes'

#Try normalization by 90 quantile
mirseq_norm = mirseq_norm/np.percentile(mirseq_norm, 90, axis=0)
#Normalize by total expression
#mirseq_norm = mirseq_norm/mirseq_norm.sum()*max(sum(mirseq_norm))
#Scale the arrays by mean and std.
mirseq_norm = MicroArray.scale(mirseq_norm)
u,s,vt=MicroArray.QaD_SVD(mirseq_norm, mirseq_meta, plotEigengenes=False)

#Filter down features to high variance ones and perform hierarchical clustering
variance = mirseq_norm.var(axis=1)
x = mirseq_norm[variance>np.percentile(variance, 95)]
clust.heatmap(x.as_matrix(), list(x.index), list(x.columns), filename='foo.png')

# <markdowncell>

# ###RNA-Seq

# <codecell>

##Remove a sample with problematic karyotype
rnaseq = rnaseq.drop(['SC12-040.420.12.19', 'SC11-008DE.134.5.18'], 1)
rnaseq_meta = rnaseq_meta.drop(['SC12-040.420.12.19', 'SC11-008DE.134.5.18'], 0)

##Remove non protein coding genes
biomart = pd.DataFrame.from_csv(syn.get(BIOMART_ANNOT_ID).path, sep='\t')
idx = biomart.index[biomart['Gene Biotype']=='protein_coding']
#Filter use first index of multindex hence i[0] remove the Ensemble id text behind '.'  
exprData = rnaseq.select(lambda i: i[0].split('.')[0] in idx)
print 'RNA-Seq was filtered from', len(rnaseq), 'features to', len(exprData), 'protein coding genes'

#Filter out genes with close to zero expression
idx = sum(exprData<=3, 1)/exprData.shape[1]>.90
print 'There were %i genes with 99%% of samples having RPKM < 3 of these' %sum(idx)
print '%i genes the mean expression was %0.2g and median %0.2g' %(sum(idx), np.mean(np.mean(exprData[idx])), 
                                                                  np.median(exprData[idx]))
exprData = exprData[~idx]
print exprData.shape

#Log Transform data and standardize
logExprData = np.log2(exprData.ix[:,:]+.00001)
logExprData = MicroArray.scale(logExprData)

#Remove control samples and look at PCs
idx = (rnaseq_meta['GroupLevel1 DifferentiationState']!='PrimaryCellCulture')
d=logExprData.ix[:,idx]
d_meta = rnaseq_meta.ix[idx,:]
u,s,vt=MicroArray.QaD_SVD(d, d_meta, plotEigengenes=False)

# <markdowncell>

# Observations: Once we drop the control samples the first through 4th PC capture EB vs SC, EB vs SC, SC/EB vs MESO, ECTO vs rest respectively.  Together these capture 85% of the information.  PC1 seems to have an outlier (SC11-008DE.134.5.18) but removing this sample shows the same trend which seems to indicate that this sample is not driving the first PC just project onto it strongly.

# <markdowncell>

# ####RNA-Seq Enrichment

# <codecell>

symbol = d.index.get_level_values(1)
idx2 = np.argsort(u[:,1])
idx3 = np.argsort(u[:,2])
idx4 = np.argsort(u[:,3])
HTML('<table><tr><td>'+
     ToppGeneEnrichement(symbol[idx2][:100], name='Correlated PC2')._repr_html_()+'</td><td>'+
     ToppGeneEnrichement(symbol[idx2][-100:], name='Anti-Corr PC2')._repr_html_()+'</td></tr><tr><td>'+
     ToppGeneEnrichement(symbol[idx3][:100], name='Correlated PC3')._repr_html_()+'</td><td>'+
     ToppGeneEnrichement(symbol[idx3][-100:], name='Anti-Corr PC3')._repr_html_()+'</td></tr><tr><td>'+
     ToppGeneEnrichement(symbol[idx4][:100], name='Correlated PC4')._repr_html_()+'</td><td>'+
     ToppGeneEnrichement(symbol[idx4][-100:], name='Anti-Corr PC4')._repr_html_()+'</td></tr></table>')

# <markdowncell>

# ###Methylation
# The methylation data has arelady been quantile normalized.  (See the provenance of the [syn2233188](https://www.synapse.org/#!Synapse:syn2233188)).  I will perform a PCA anyway to look for issues with the data by first extracting the top 50th percentile variable probes.

# <codecell>

#The methylation data has arelady been quantile normalized will perform 
#scaling before performing PCA and reduce the dataset (for speed) by
var = methyl.var(axis=1)
d = methyl[var>np.percentile(var, 50)]
d = MicroArray.scale(d)
u,s,vt=MicroArray.QaD_SVD(d, methyl_meta, plotEigengenes=False)

# <markdowncell>

# **Observations:** The PCs correspond well with differentiation state.  Scaling by features (not shown) shows that the 4th pc corresponds to Gender.

# <markdowncell>

# #Save normalized samples to Synapse

# <codecell>

code = synapseclient.File('C4_qc_filtering_normalization.ipynb', parentId='syn2246673')
code = syn.store(code)

logExprData.to_csv('rnaseq_norm.tsv', sep='\t')
file = synapseclient.File('rnaseq_norm.tsv', parentId='syn1773109', name='Normalized mRNACalls')
rnaID = syn.store(file, used=[BIOMART_ANNOT_ID, EXPR_ID, EXPR_META_ID], 
                  executed=[code],
                  activityName='QC and Filtering')

#mirseq_norm.to_csv('mirseq_norm.tsv', sep='\t')
#mirID = syn.store(synapseclient.File('mirseq_norm.tsv', parentId='syn1773109', name='Normalized miRNACalls'))

# <markdowncell>

# #Follow-up analysis
# is available in the notebook: C4-differential_analysis.ipynb

