# -*- coding: utf-8 -*-
"""
Created on Mon May 23 22:07:30 2022

@author: ftead
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import math
import re
from matplotlib_venn import venn2
import gseapy as gp

#%%

#Took the data that already has all the differentially expressed genes.
#The data already has DE genes and segnificant padj genes.

data = pd.read_csv('WT_NDvsKO_ND_DE.csv', sep = '\s+', engine = 'python')
#Filtering only protein coding genes. 
data = data.loc[data['GENEBIOTYPE'] == 'protein_coding']
#data['AllCAPS'] = data['SYMBOL'].str.upper()
#%%

#Enrichment analysis librarys names
libs = gp.get_library_name(organism="Mouse")
gene_sets = [
             'ChEA_2016','ENCODE_TF_ChIP-seq_2015','ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X','GO_Biological_Process_2021',
             'GO_Cellular_Component_2021','GO_Molecular_Function_2021','GTEx_Aging_Signatures_2021','GWAS_Catalog_2019',
             'KEGG_2019_Mouse','KEGG_2021_Human','Reactome_2016'
             ]
#%%
#differentiating to upregulated and down regulated genes
data_up = data.loc[data['log2FoldChange'] > 0]
data_down = data.loc[data['log2FoldChange'] < 0]
#%%

#Performing ennrichment analysis for upregulated, downregulated and original database.
enr_all = gp.enrichr(gene_list = list(data['SYMBOL'].unique()), description='test_name', gene_sets= gene_sets, cutoff = 0.5)
enr_up = gp.enrichr(gene_list = list(data_up['SYMBOL'].unique()), description='test_name', gene_sets= gene_sets, cutoff = 0.5)
enr_down = gp.enrichr(gene_list = list(data_down['SYMBOL'].unique()), description='test_name', gene_sets= gene_sets, cutoff = 0.5)

results_all = enr_all.results
results_up = enr_up.results
results_down = enr_down.results

results_all.to_csv('enricher_results_all.csv')
results_up.to_csv('enricher_results_up.csv')
results_down.to_csv('enricher_results_down.csv')

#%%

#Importing atac-seq data
atac = pd.read_excel('dif_acces_genes_narrowpeak_atacseq_fc_padj.xlsx', index_col = 0)

#%%

#Enrichment analysis of atac seq data

enr_atac = gp.enrichr(gene_list = list(atac['SYMBOL'].unique()), description='test_name', gene_sets= gene_sets, cutoff = 0.5)
results_atac = enr_atac.results
results_atac.to_csv('enricher_results_atac.csv')

#%%

#Filtering all the segnificantly enriched categories
results_all_seg = results_all.loc[results_all['Adjusted P-value']<0.05]
results_up_seg = results_up.loc[results_up['Adjusted P-value']<0.05]
results_down_seg = results_down.loc[results_down['Adjusted P-value']<0.05]
results_atac_seg = results_atac.loc[results_atac['Adjusted P-value']<0.05]

#%%

#preforming hypergeometric test of the cross between the RNA-seq upregulated genes and atac-seq genes

#Choosing data onlt for ChIPX data

results_all_seg_ChIPX = results_all_seg.loc[results_all_seg['Gene_set'] == 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X']
results_up_seg_ChIPX = results_up_seg.loc[results_up_seg['Gene_set'] == 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X']
results_down_seg_ChIPX = results_down_seg.loc[results_down_seg['Gene_set'] == 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X']
results_atac_seg_ChIPX = results_atac_seg.loc[results_atac_seg['Gene_set'] == 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X']

pvval_hgeom = stats.hypergeom.pmf(len(results_atac_seg_ChIPX['Term'][results_atac_seg_ChIPX['Term'].isin(results_up_seg_ChIPX['Term'])]),104,len(results_atac_seg_ChIPX['Term']), len(results_up_seg_ChIPX['Term']))

#%%

#producing venn (Figure 1 B)

x = len(results_atac_seg_ChIPX['Term'][results_atac_seg_ChIPX['Term'].isin(results_up_seg_ChIPX['Term'])])
a = len(results_atac_seg_ChIPX['Term'])
r = len(results_up_seg_ChIPX['Term'])

plt.figure(figsize=(4,5), dpi = 600)
venn2(subsets = (a-x, r-x,x),set_labels = None,set_colors=('pink', 'blue'), alpha = 0.4)
plt.savefig('ATAC_RNA_venn.png',bbox_inches='tight',dpi =300)

#%%

#orgnizing data for enrichment plot
go_terms = ['GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021']

results_all_seg_GO = results_all_seg[results_all_seg.Gene_set.isin(go_terms)]
results_up_seg_GO = results_up_seg[results_up_seg.Gene_set.isin(go_terms)]
results_down_seg_GO = results_down_seg[results_down_seg.Gene_set.isin(go_terms)]

#%%
df = pd.DataFrame()

for i in go_terms:
    df = pd.concat([df,results_up_seg_GO.loc[results_up_seg_GO['Gene_set'] == i].iloc[0:5]])
    
df['-log(Adjusted P-value)'] =np.log10(df['Adjusted P-value'])*-1
#%%

#plotting the figure (Figure 1 A)

fig, ax = plt.subplots()
plt.figure(figsize=(1,8), dpi = 300)
sns.set_theme()
sns.set(font_scale = 1.4)
sns.set_style("whitegrid")
ax = sns.relplot(data = df, x = '-log(Adjusted P-value)', y = 'Term', hue = 'Gene_set',size ='Odds Ratio',sizes=(50,300), palette = 'mako',edgecolor="black")
plt.xlim(0)
plt.savefig('go_results_up_genes.png',dpi = 300,bbox_inches='tight',transparent=True)


#%%

#Generating correlation matrix with SIRT6 and REST, SMAD4 and SUZ12 from allen brain atlas brains.

#Impoering and organizing the DATA
def Expression_to_Columns(Columns_file, Expression_file, Gene):
    df_Columns = pd.read_csv(Columns_file)
    df_Columns = df_Columns.replace(-1,np.NaN).dropna(axis = 0)
    df_Expression = pd.read_csv(Expression_file, index_col = 0, header = None)
    df_Expression = df_Expression.transpose().reset_index().drop(columns = 'index').dropna()
    df_Expression = df_Expression.rename(columns={ df_Expression.columns[0]: 'Expression' + ' ' + Gene})
    df_col_exp = pd.concat([df_Columns,df_Expression], axis = 1)
    return df_col_exp

joined_REST = Expression_to_Columns('REST_Columns.csv', 'REST_Expression.csv', 'REST')
joined_SIRT6 = Expression_to_Columns('SIRT6_Columns.csv', 'SIRT6_Expression.csv', 'SIRT6')
joined_SUZ12 = Expression_to_Columns('SUZ12_Columns.csv', 'SUZ12_Expression.csv', 'SUZ12')
joined_SMAD4 = Expression_to_Columns('SMAD4_Columns.csv', 'SMAD4_Expression.csv', 'SMAD4')

joined = pd.concat([joined_REST, joined_SIRT6['Expression SIRT6'], joined_SUZ12['Expression SUZ12'], joined_SMAD4['Expression SMAD4']], axis = 1)

#%%
#Correlation


corr_matrix = pd.DataFrame()
for i in list(joined['donor_age'].unique()):
    j = joined.loc[joined['donor_age'] == i]
    s6_rest = round(j['Expression SIRT6'].corr(j['Expression REST']),3)
    s6_suz12 = round(j['Expression SIRT6'].corr(j['Expression SUZ12']),3)
    s6_smad4 = round(j['Expression SIRT6'].corr(j['Expression SMAD4']),3)
    corr = pd.DataFrame({'Age':[i], 'REST': [s6_rest], 'SUZ12': [s6_suz12], 'SMAD4': [s6_smad4]})
    corr_matrix = pd.concat([corr_matrix,corr])

corr_matrix = corr_matrix.set_index(corr_matrix.columns[0])

#%%
#plotting correlation matrix (Figure 1 C)
fig, g = plt.subplots()
plt.figure(figsize=(1,1), dpi = 300)
sns.set_theme()
sns.set(font_scale = 2)
g = sns.clustermap(corr_matrix, vmin = -1,vmax = 1,cmap="coolwarm")
plt.savefig('correlation_heatmap_tfs_vs_sirt6.png', dpi = 300, bbox_inches = 'tight', transparent = True)


#%%

#Charecteraizing REST targets
#Reading peak annotation resolts from GSM803365 H1 REST ChIP-seq Broadpeak files, generated in R
rep1 = pd.read_csv('annotated peaks REST ChIP-seq rep 1 GSM803365.csv',
                   sep = ',', engine = 'python', index_col = False).dropna().drop_duplicates(subset=['SYMBOL'])
rep2 = pd.read_csv('annotated peaks REST ChIP-seq rep 2 GSM803365.csv',
                   sep = ',', engine = 'python',index_col = False).dropna().drop_duplicates(subset=['SYMBOL'])

#%%

rest_unique_targets_in_h1 = pd.Series(list(rep1['SYMBOL'][rep1['SYMBOL'].isin(rep2['SYMBOL'])]))
data_up['Upper'] = data_up['SYMBOL'].str.upper()

ko_up = data_up['Upper']

inter = pd.Series(list(set(rest_unique_targets_in_h1) & set(ko_up)))
#%%
pvval_hgeom_restT_upreg = stats.hypergeom.pmf(len(inter),20203,len(ko_up), len(rest_unique_targets_in_h1))


plt.figure(figsize=(4,5), dpi = 600)
font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 18}
plt.rc('font', **font)
out = venn2(subsets = (len(ko_up) - len(inter), len(rest_unique_targets_in_h1) - len(inter),len(inter)),set_labels = None,set_colors=('pink', 'blue'), alpha = 0.4)
plt.savefig('upregulated in ko vs REST targets venn.png',bbox_inches='tight',dpi =300)
#%%
#Enrichment analysis for the intersecting genes
inter_enr = gp.enrichr(gene_list = inter, description='test_name', gene_sets= gene_sets, cutoff = 0.5)
results_inter = inter_enr.results
#%%
results_inter.to_csv('upregulated REST targents enrichmnet.csv')
#%%
results_inter_seg = results_inter.loc[results_inter['Adjusted P-value'] < 0.05]
results_inter_seg_GO = results_inter_seg[results_inter_seg.Gene_set.isin(go_terms)]

df2 = pd.DataFrame()

for i in go_terms:
    df2 = pd.concat([df2,results_inter_seg_GO.loc[results_inter_seg_GO['Gene_set'] == i].iloc[0:5]])
    
df2['-log(Adjusted P-value)'] =np.log10(df2['Adjusted P-value'])*-1

#%%

fig, ax = plt.subplots()
plt.figure(figsize=(1,8), dpi = 300)
sns.set_theme()
sns.set(font_scale = 1.4)
sns.set_style("whitegrid")
ax = sns.relplot(data = df2, x = '-log(Adjusted P-value)', y = 'Term', hue = 'Gene_set',size ='Odds Ratio',sizes=(50,300), palette = 'mako',edgecolor="black")
plt.xlim(0)
plt.savefig('go_results_up_genes_rest_targets.png',dpi = 300,bbox_inches='tight',transparent=True)
