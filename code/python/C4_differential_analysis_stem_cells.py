# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import synapseclient
import numpy as np
import pandas as pd
import scipy.stats as stats
#Regression 
import statsmodels.formula.api as sm
import statsmodels.stats.multitest as smm
from statsmodels.regression.mixed_linear_model import MixedLM
#Analysis tools
from collections import Counter as count
import hierarchical_clustering as clust
import MicroArray
from toppGenePost import ToppGeneEnrichement
#Other packages
import IPython.display as display
import string

EXPR_NORM_ID = 'syn2701943'
EXPR_META_ID = 'syn2278178'
METH_ID = 'syn2233188'
METH_META_ID = 'syn2677043'
MIRSEQ_NORM_ID = 'syn2701942'
MIRSEQ_META_ID = 'syn2278179'
pd.set_option('max_columns', 200)
syn = synapseclient.login()

# <markdowncell>

# #Download and annotate data
# 1. Downloads expression/miRNA-Seq/methylation data from Synapse
# 2. Gets annotations from metadata in Synapse
# 3. Make sure the metadata matches up with samples

# <codecell>

mirseq = pd.read_csv(syn.get(MIRSEQ_NORM_ID).path, sep='\t', index_col=0)
mirseq_meta = pd.read_csv(syn.get(MIRSEQ_META_ID).path, sep='\t',index_col=0)
mirseq_meta = mirseq_meta.ix[mirseq.columns,:]
print 'Mirseq Data Size:', mirseq.shape, mirseq_meta.shape
rnaseq = pd.read_csv(syn.get(EXPR_NORM_ID).path, sep='\t', index_col=[0,1,2])
rnaseq_meta = pd.read_csv(syn.get(EXPR_META_ID).path, sep='\t',index_col=0)
rnaseq_meta = rnaseq_meta.ix[rnaseq.columns,:]
print 'RNA-Seq Data Size:', rnaseq.shape, rnaseq_meta.shape
methyl = pd.read_csv(syn.get(METH_ID).path, sep='\t', index_col=0)
methyl_meta = pd.read_csv(syn.get(METH_META_ID).path, sep='\t',index_col=0)
methyl_meta = methyl_meta.ix[methyl.columns,:]
print 'Methylation Data Size:', methyl.shape, methyl_meta.shape

# <markdowncell>

# ###Prune and rename metadata

# <codecell>

toRemove = ['Related C4 ID comment', 'Sage verified RNA-Seq', 'Sage verified']
mirseq_meta.drop(toRemove,1, inplace=True)
rnaseq_meta.drop(toRemove,1, inplace=True)
rnaseq_meta.drop(['Bad Lines'],1, inplace=True)
methyl_meta.drop(toRemove,1, inplace=True)
mirseq_meta.rename(columns={'High Confidence Donor ID (HCDID)':'Donor id'}, inplace=True)
rnaseq_meta.rename(columns={'High Confidence Donor ID (HCDID)':'Donor id'}, inplace=True)
methyl_meta.rename(columns={'High Confidence Donor ID (HCDID)':'Donor id'}, inplace=True)
methyl_meta['Array'] = methyl_meta['Array'].astype('str')

# <markdowncell>

# ###Filter down to the stem cell only samples

# <codecell>

mirseq = mirseq.ix[:,mirseq_meta.CellDiffState=='SC']
mirseq_meta = mirseq_meta.ix[mirseq_meta.CellDiffState=='SC',:]
print 'Mirseq Data Size:', mirseq.shape, mirseq_meta.shape
rnaseq = rnaseq.ix[:, rnaseq_meta.CellType=='SC']
rnaseq_meta = rnaseq_meta.ix[rnaseq_meta.CellType=='SC',:]
print 'RNA-seq Data Size:', rnaseq.shape, rnaseq_meta.shape
methyl = methyl.ix[:,methyl_meta['Diffname short']=='SC']
methyl_meta = methyl_meta.ix[methyl_meta['Diffname short']=='SC', :]
print 'Methylation Data Size:', methyl.shape, methyl_meta.shape

# <markdowncell>

# #Make PCA plots

# <markdowncell>

# ##RNA-Seq

# <codecell>

u,s,vt = MicroArray.QaD_SVD(MicroArray.scale(rnaseq.T).T, rnaseq_meta.drop('CellLine', 1), plotEigengenes=False)

# <markdowncell>

# ##Micro RNA-Seq

# <codecell>

u,s,vt = MicroArray.QaD_SVD(MicroArray.scale(mirseq.T).T, mirseq_meta, plotEigengenes=False)

# <markdowncell>

# ##Methylation Arrays

# <codecell>

u,s,vt = MicroArray.QaD_SVD(MicroArray.scale(methyl.T).T, methyl_meta, plotEigengenes=False)

# <markdowncell>

# #Identify features correlated with reprogramming and cell of origin
# 
# As we have multiple samples derived from the same individual we will used a mixed Effect Linear model where we treat donor and gender as radom effects.  In effect we are fitting the following model for each gene/feature $i$.
# $$Y_{i} = X\cdot\beta + Z\cdot\gamma + epsilon$$
# 
# where $X$ represent for example the  reprogramming vector or cell of origin an $Z$ is the matrix of random effects components.
# 

# <codecell>

def lmemodel(data, metadata, 
             fixedEffects = ['Tissue of Origin'],
             randomEffects=['High Confidence Donor ID (HCDID)']):
    """Performs a mixed effect linear model"""
    df = metadata[fixedEffects].copy()
    df = pd.concat([df, metadata[randomEffects]], axis=1 )
    #Change the parameters to be compatible with patsy formulas
    fixedEffects = [c.translate(string.maketrans(' ()', '___')) for c in fixedEffects]
    randomEffects = [c.translate(string.maketrans(' ()', '___')) for c in randomEffects]
    df.columns = [c.translate(string.maketrans(' ()', '___')) for c in df.columns]

    model_string = 'gene ~ '+' + '.join(fixedEffects)
    results = []
    for i in range(data.shape[0]):
        #Add the dependent variable to the dataframe
        df['gene'] = data.irow(i)
        #################
        df['High_Confidence_Donor_ID__HCDID_'] = stats.binom.rvs(1, .4, size=69)
        print df.shape, model_string
        return df
        df = df.dropna()
        print df.shape
        #################
        #compute new model
        mod = MixedLM.from_formula(model_string, df, groups = df[randomEffects])
        return mod
        #df.boxplot(by=fixedEffects)
        #mod = sm.ols(model_fit, df)
        #results.append(mod.fit())
    return results



#pd.crosstab(rnaseq_meta['Tissue of Origin'], rnaseq_meta['High Confidence Donor ID (HCDID)'])
res = lmemodel(rnaseq.ix[:1,:], rnaseq_meta)
#Find the significant genes based on f-test
#pvals = np.asarray([r.f_pvalue for r in res])
#print sum(pvals*len(pvals)<0.01), pvals

#Plot some of the genes
#pd.concat( [mirseq.ix[:1,:].T, mirseq_meta[['Tissue of Origin', 
#                                            'Gender', 
#                                            'High Confidence Donor ID (HCDID)' ]]], axis=1).boxplot(by=['Tissue of Origin', 
#                                                                                                        'Gender'], rot=90)

# <codecell>

#res.fit()
pd.crosstab(res.Tissue_of_Origin, res.High_Confidence_Donor_ID__HCDID_)

# <codecell>

def residualmodel(data, metadata, 
          parameters=['Gender', 'High Confidence Donor ID (HCDID)', 
                      'Tissue of Origin']):
    """Performs a lineal mode on all parameters except last one then uses the residuals of this model
    to fit a model on the last parameters"""
    df = metadata[parameters].copy()
    #Change the parameters to be compatible with patsy formulas
    df.columns = [c.translate(string.maketrans(' ()', '___')) for c in df.columns]
    parameters = list(df.columns)
    norm_params = parameters[:-1]
    fit_params = parameters[-1]
    model_norm = 'gene ~ '+' + '.join(norm_params)
    model_fit = 'gene ~ '+ fit_params
    results = []
    for i in range(data.shape[0]):
        df['gene'] = data.irow(i)
        #Fit model for nuiscance parameters
        mod = sm.ols(model_norm, df)
        #Find residuals and compute new model
        df['gene'] = mod.fit().resid
        df.boxplot(by=parameters[-1])
        mod = sm.ols(model_fit, df)
        results.append(mod.fit())
    return results

def model(data, metadata, parameter='Gender'):
    """Performs a lineal mode on all parameters except last one then uses the residuals of this model
    to fit a model on the last parameters"""
    df = metadata[[parameter]].copy()
    #Change the parameters to be compatible with patsy formulas
    df.columns = [parameter.translate(string.maketrans(' ()', '___'))]
    parameter = df.columns[0]
    model_fit = 'gene ~ '+ parameter
    results = []
    for i in range(data.shape[0]):
        df['gene'] = data.irow(i)
        #Fit model for nuiscance parameters
        mod = sm.ols(model_fit, df)
        results.append(mod.fit())
    return results

resOrigin = model(mirseq.ix[:,:], mirseq_meta, parameter='Tissue of Origin')
resGender = model(mirseq.ix[:,:], mirseq_meta, parameter='Gender')
resReprog = model(mirseq.ix[:,:], mirseq_meta, parameter='Reprogramming Vector Type')
resCellType = model(mirseq.ix[:,:], mirseq_meta, parameter='Cell Line Type')
pvals = pd.DataFrame({'Tissue_of_Origin':[r.f_pvalue for r in resOrigin],
                      'Gender' :[r.f_pvalue for r in resGender],
                      'Reprogramming_Vecor': [r.f_pvalue for r in resReprog],
                      'Cell Line Type': [r.f_pvalue for r in resCellType]},
                     index = mirseq.index)

# <codecell>

pvals.sort(pvals.columns[3])
#q=pd.concat([mirseq_meta, mirseq.ix[pvals[4:6].index,:].T], axis=1)
#q.boxplot(by='Reprogramming Vector Type')
#mirseq_meta[mirseq_meta['Reprogramming Vector Type']=='lentivirus']
#pd.crosstab(mirseq_meta['Reprogramming Vector Type'], mirseq_meta['High Confidence Donor ID (HCDID)'])

# <codecell>

mirseq_meta

# <markdowncell>

# ##Cell of origin and Induction gene analysis

# <codecell>

def runModels(data):
    #Put the data into temporary datastructure
    df = data.T
    #I have to temporarily rename the genes to work with Patsy formulas
    df.columns = ['gene_%i' %i for i in range(len(df.columns))] 
    df['origcell'] = metadata['origcell']
    df['inductiongenes'] = metadata['inductiongenes']

    #Exclude HESCs
    #idx = metadata.linetype !='hESC'
    #df = df.ix[idx,:]
    
    def model(geneName):
        mod1 = sm.ols('%s ~ origcell' %geneName, df)
        mod2 = sm.ols('%s ~ inductiongenes' %geneName, df)
        return mod1.fit().f_pvalue, mod2.fit().f_pvalue
    
    #from multiprocessing import Pool
    #pool = Pool()
    pvals = map(model, df.columns[:-2])
    #pool.join()
    pvals = np.asarray(pvals)
    
    models = pd.DataFrame(symbol.copy())
    models['origcell']=pvals[:,0]
    models['inductiongenes']=pvals[:,1]
    return models

normPvals = runModels(logExprDataNorm)
#origPvals = runModels(logExprData)

# <rawcell>

# ##Store the pvalues to synapse
# models = pd.DataFrame(symbol.copy())
# models['coi_fdr_corrected_p (normalized data)'] = smm.multipletests(normPvals.origcell, alpha=0.05, method='fdr_bh')[1]
# models['genes_fdr_corrected_p (normalized data)'] = smm.multipletests(normPvals.inductiongenes, alpha=0.05, method='fdr_bh')[1]
# 
# models['coi_fdr_corrected_p (uncorrected data)'] = smm.multipletests(origPvals.origcell, alpha=0.05, method='fdr_bh')[1]
# models['genes_fdr_corrected_p (uncorrected data)'] = smm.multipletests(origPvals.inductiongenes, alpha=0.05, method='fdr_bh')[1]
# 
# models.to_csv('model_corrected_pvalues_for_shane.csv')
# 
# code = syn.store(synapseclient.File('C4_qc_and_linear_models.ipynb', parent='syn1774100'))
# ent = synapseclient.File('model_corrected_pvalues_for_shane.csv', parent='syn2332184')
# ent = syn.store(ent, used=[EXPR_ID, METADATA_ID], executed=code)
# syn.onweb(ent)   

# <markdowncell>

# ###Compare these models to the genelists found by Nathan

# <codecell>

import matplotlib_venn

#Determine significant genes using Benjamini-Hochberg FDR correction
def filterPvals2GeneSymbols(pvals):
    pval_pass = smm.multipletests(pvals, alpha=0.05, method='fdr_bh')[0]
    return symbol[pval_pass]

#coi_signficant_logExprData= filterPvals2GeneSymbols(origPvals.origcell)
coi_signficant_logExprDataNorm= filterPvals2GeneSymbols(normPvals.origcell)
nathan_coi=pd.DataFrame.from_csv('COI-GeneExpression-signatures.csv')
matplotlib_venn.venn2([set(coi_signficant_logExprDataNorm),  #set(coi_signficant_logExprData), 
                       set(nathan_coi.index)], 
                       set_labels=['normalized', 'AltAnalyze']) #'unormalized', 
pylab.title('Significant Gene counts for COI')

display.display_html(ToppGeneEnrichement(coi_signficant_logExprDataNorm, name='COI'))
print set(nathan_coi.index).intersection(set(coi_signficant_logExprDataNorm))


# <codecell>

#inductiongenes_signficant_logExprData= filterPvals2GeneSymbols(origPvals.inductiongenes)
inductiongenes_signficant_logExprDataNorm= filterPvals2GeneSymbols(normPvals.inductiongenes)
nathan_inductiongenes=pd.DataFrame.from_csv('GeneCombinations-GeneExpression-signatures.csv')
matplotlib_venn.venn2([set(inductiongenes_signficant_logExprDataNorm), #set(inductiongenes_signficant_logExprData), 
                       set(nathan_inductiongenes.index)], 
                       set_labels=['normalized', 'AltAnalyze'])
pylab.title('Significant Gene counts for Induction Genes')

# <markdowncell>

# ####Explore some of the specific different genes

# <codecell>


# <markdowncell>

# **Ah** This seems to be due to HESC vs iPSC - should HESC be excluded for COI analysis?

