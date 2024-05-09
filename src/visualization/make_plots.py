import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib_venn as venn

# LOAD DATA
obs = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/observations_HL60.csv')
kinases = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/r_things/kpas_data/kinases_final.csv')
exp_in = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/expressed_in.csv')
pk_rel = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pert_kin_rel.csv')
ks_rel = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/enz_sub_omnipath.csv')
ksea_df = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/ksea_HL60.csv')
pdts = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/pdts_HL60.csv')
lints = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/lints_HL60.csv')
lints_nscorr = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/lints_nscorr_HL60.csv')
dpdts = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/dirts_pdt_HL60.csv')
dirts = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/dirts_HL60.csv')
dirts_nscorr = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/dirts_nscorr_HL60.csv')


# OBSERVATIONS
tprot_tup = set(obs['tprot'].values.tolist())
kinases_tup = set(kinases['kinase'].values.tolist())
exp_tup = set(exp_in.loc[exp_in['expressed_in'] == 'HL60', 'kinase_name'].values.tolist())
venn.venn3_unweighted([tprot_tup, kinases_tup, exp_tup], set_labels=('tprot', 'kinase', 'exp_in'))


# KNOWN INHIBITORS
# perturbagen/kinase coverage
# histogram: perturbagens per kinase / kinases per perturbagen
p_count = pk_rel['perturbagen'].value_counts()
k_count = pk_rel['kinase'].value_counts()
# plot
plt.hist(p_count, bins=200)
plt.xlabel('nr. of kinases')
plt.ylabel('perturbagens')
plt.hist(k_count)
plt.xlabel('nr. of perturbagens')
plt.ylabel('kinases')

# heatmap
pk_hm = pk_rel[['kinase', 'perturbagen']].copy()
pk_hm['score'] = pk_rel[['discoverx', 'kuster']].mean(axis=1)
# data is converted into df, where rows and cols correspond to phosphosite and perturbagens, respectively.
pk_hm = pk_hm.pivot(index='perturbagen', columns='kinase', values='score')
pk_hm = pk_hm.replace(np.nan, 1)
# plot
plt.pcolor(pk_hm)
plt.yticks(np.arange(0.5, len(pk_hm.index), 1), pk_hm.index)
plt.xticks(np.arange(0.5, len(pk_hm.columns), 1), pk_hm.columns)
plt.xticks(rotation=90)
plt.colorbar()

# kinase coverage
k_pert = set(pk_rel['kinase'])
venn.venn3_unweighted([tprot_tup, kinases_tup, k_pert], set_labels=('tprot', 'kinase', 'k_pert'))
venn.venn3_unweighted([exp_tup, kinases_tup, k_pert], set_labels=('exp_in', 'kinase', 'k_pert'))

# KSEA
# score distributions
# histogram
ksea_df_tc = ksea_df.groupby('kpa')['tc'].mean()
ksea_df_tc = ksea_df_tc[ksea_df_tc <= 50]
#ksea_df_tc = ksea_df.loc[ksea_df['tc'] <= 50]
n, bins, patches = plt.hist(x=ksea_df_tc, bins=50, color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('# of targets')
plt.ylabel('Frequency')
maxfreq = n.max()

# Barplot
ksea_df['status'].value_counts().plot(kind='bar')


# PDTs, LINTs and DIRTs
# tailor known_targets data to cell line
ks_rel = ks_rel[ks_rel['kpa'].isin(pk_rel['kinase'])]
ks_rel = ks_rel[ks_rel['pst'].isin(obs['substrate'])]
ks_rel = ks_rel[['kpa', 'pst']]

# venn diagram ks links
ks_tup = set(tuple(i) for i in ks_rel.values.tolist())

pdts_tup = set(tuple(j) for j in pdts[['kinase', 'substrate']].values.tolist())
lints_nscorr_tup = set(tuple(k) for k in lints_nscorr[['kinase', 'substrate']].values.tolist())
lints_tup = set(tuple(m) for m in lints[['kinase', 'substrate']].values.tolist())

dpdts_tup = set(tuple(n) for n in dpdts[['K', 'PsT']].values.tolist())
dirts_nscorr_tup = set(tuple(o) for o in dirts_nscorr[['K', 'PsT']].values.tolist())
dirts_tup = set(tuple(p) for p in dirts[['K', 'PsT']].values.tolist())

venn.venn3_unweighted([ks_tup, pdts_tup, lints_nscorr_tup], set_labels=('ks', 'pdt', 'lint_nscorr'))
venn.venn3_unweighted([ks_tup, lints_tup, lints_nscorr_tup], set_labels=('ks', 'lint', 'lint_nscorr'))
venn.venn3_unweighted([ks_tup, dpdts_tup, dirts_nscorr_tup], set_labels=('ks', 'dpdt', 'dirt_nscorr'))
venn.venn3_unweighted([ks_tup, dirts_tup, dirts_nscorr_tup], set_labels=('ks', 'dirt', 'dirt_nscorr'))

dpdts_MCF7 = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/dirts_pdt_MCF7.csv')
dpdts_NTERA2 = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/dirts_pdt_NTERA2.csv')
dirts_nscorr_MCF7 = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/dirts_nscorr_MCF7.csv')
dirts_nscorr_NTERA2 = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdt_approach/dirts_nscorr_NTERA2.csv')

dpdts_MCF7_tup = set(tuple(n) for n in dpdts_MCF7[['K', 'PsT']].values.tolist())
dpdts_NTERA2_tup = set(tuple(n) for n in dpdts_NTERA2[['K', 'PsT']].values.tolist())
dirts_MCF7_nscorr_tup = set(tuple(o) for o in dirts_nscorr_MCF7[['K', 'PsT']].values.tolist())
dirts_NTERA2_nscorr_tup = set(tuple(o) for o in dirts_nscorr_NTERA2[['K', 'PsT']].values.tolist())

venn.venn3_unweighted([dpdts_tup, dpdts_MCF7_tup, dpdts_NTERA2_tup], set_labels=('dpdt_HL60', 'dpdt_MCF7', 'dpdt_NTERA2'))
venn.venn3_unweighted([dirts_nscorr_tup, dirts_MCF7_nscorr_tup, dirts_NTERA2_nscorr_tup], set_labels=('dirt_ns_HL60', 'dirt_ns_MCF7', 'dirt_ns_NTERA2'))
