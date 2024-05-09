import pandas as pd
import numpy as np
from pyScripts.lm_metrics import explain_psts

# INPUT
pert_kin = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pert_kin_rel.csv')
# EDGES
enz_tprot = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/enz_tprot_strong_HL60.csv')
enz_sub = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/enz_sub_strong_HL60.csv')
# OUTPUT
obs = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/observations_HL60.csv')
# UNIVERSE
enz_sub_univ = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/enz_sub_omnipath.csv')


# perturbagens IN
perturbagens = pert_kin['perturbagen'].drop_duplicates().tolist()
pert_in = 'Torin'  # select one

# subsetting data
# perturbed (affected) psts OUT
obs_pert = obs.loc[(abs(obs['fold_change']) > 0) &
                   (obs['p_value'] <= 0.05) &
                   (obs['perturbagen'] == pert_in)].reset_index(drop=True)
# edges
enz_sub = enz_sub.loc[enz_sub['Pert'] == pert_in].reset_index(drop=True)
enz_tprot = enz_tprot.loc[enz_tprot['Pert'] == pert_in].reset_index(drop=True)
# psts KNOWN in enzyme-substrate relationship(s) universe
psts_known = list(set(enz_sub_univ['pst']))


# EXPLAIN
# edges_df
pert_kin_edges = np.array(pert_kin.loc[pert_kin['perturbagen'] == pert_in, ['perturbagen', 'kinase']])
enz_tprot_edges = np.array(enz_tprot.loc[enz_tprot['Pert'] == pert_in, ['Kpa', 'Tprot']])
enz_sub_edges = np.array(enz_sub.loc[enz_sub['Pert'] == pert_in, ['Kpa', 'Pst']])
edges_df = pd.DataFrame(np.concatenate((pert_kin_edges, enz_tprot_edges, enz_sub_edges)),
                        columns=['from', 'to'])

# psts affected by perturbagen
psts_pert = obs_pert.loc[obs_pert['perturbagen'] == pert_in, 'substrate'].tolist()

# determine explainable phosphosites and explainability
psts_expl, expl = explain_psts(pert_in, edges_df, psts_pert, psts_known)

expl_all_psts = expl[4]
expl_known_psts = expl[5]


# SUMMARY
import matplotlib_venn as venn

tprot_explained = set(enz_sub.loc[enz_sub['Pst'].isin(psts_expl), 'Tprot'])
tprot_pert = set(obs_pert.loc[obs_pert['perturbagen'] == pert_in, 'tprot'].tolist())
tprot_known = set(tprot_pert) & set(enz_sub_univ['tprot'])

len(tprot_explained) / len(tprot_pert)
len(tprot_explained) / len(tprot_known)

venn.venn3_unweighted([set(psts_pert), set(enz_sub_univ['pst']), set(psts_expl)], set_labels=('psts_pert', 'enz_sub', 'psts_explained'))
venn.venn3_unweighted([tprot_pert, set(enz_sub_univ['tprot']), tprot_explained], set_labels=('tprot_pert', 'enz_tprot', 'tprot_explained'))
