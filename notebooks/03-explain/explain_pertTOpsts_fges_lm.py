import pandas as pd
import numpy as np
from pyScripts.lm_metrics import explain_psts

# INPUT
pert_kin = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pert_kin_rel.csv')
# EDGES
enz_tprot = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/results_lm/enz_tprot_fges_30_bs10_AND_lm_MCF7.csv')
enz_sub = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/results_lm/enz_sub_fges_30_bs10_AND_lm_MCF7.csv')
# OUTPUT
obs = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/ctamdb_data/ctamdb_dpoa_MCF7.tsv', sep='\t')
# UNIVERSE
enz_sub_univ = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/fges_data/ctamdb_data/results/fges_enzsub30_bs10_MCF7.txt', sep='\t')


# perturbagens IN
perturbagens = pert_kin['perturbagen'].drop_duplicates().tolist()

# subsetting data
# perturbed (affected) psts OUT
obs_pert = obs.loc[(abs(obs['fc']) > 0) &
                   (obs['sid_score'] <= 0.05) &
                   (obs['perturbagen'].isin(perturbagens))].reset_index(drop=True)
# edges
enz_sub = enz_sub.loc[enz_sub['Pert'].isin(perturbagens)].reset_index(drop=True)
enz_tprot = enz_tprot.loc[enz_tprot['Pert'].isin(perturbagens)].reset_index(drop=True)
# psts KNOWN in enzyme-substrate relationship(s) universe
psts_known = list(set(enz_sub_univ['pst']))
# psts_known = list(set(obs.loc[obs['tprot'].isin(enz_sub_univ['tprot']), 'substrate']))  # for omnipath_tprot


# EXPLAIN
explainability = []
psts_explained = []

for pert_in in perturbagens:

    # edges_df
    pert_kin_edges = np.array(pert_kin.loc[pert_kin['perturbagen'] == pert_in, ['perturbagen', 'kinase']])
    enz_tprot_edges = np.array(enz_tprot.loc[enz_tprot['Pert'] == pert_in, ['Kpa', 'Tprot']])
    enz_sub_edges = np.array(enz_sub.loc[enz_sub['Pert'] == pert_in, ['Kpa', 'Pst']])
    edges_df = pd.DataFrame(np.concatenate((pert_kin_edges, enz_tprot_edges, enz_sub_edges)),
                            columns=['from', 'to'])

    # psts affected by perturbagen
    psts_pert = obs_pert.loc[obs_pert['perturbagen'] == pert_in, 'pst'].tolist()

    # determine explainable phosphosites and explainability
    psts_expl, expl = explain_psts(pert_in, edges_df, psts_pert, psts_known)

    # add results to dataframe
    explainability.append(expl)
    psts_explained = list(np.unique(psts_explained + psts_expl))

explainability = pd.DataFrame(explainability, columns=['pert_in', 'psts_explained', 'psts_pert', 'psts_known', 'expl_all_psts', 'expl_known_psts'])


# export to csv with:
explainability.to_csv('/Users/magdalena/Desktop/bezzlab/logic_things/fges_data/ctamdb_data/results/expl_fges_30_bs10_AND_lm_MCF7.csv', index=False)


# SUMMARY (REDO IN R)
import statistics
import matplotlib_venn as venn

all_psts_explained = set(psts_explained)
all_psts_pert = set(obs_pert['substrate'])
all_psts_known = set(all_psts_pert) & set(psts_known)

len(all_psts_explained) / len(all_psts_pert)
len(all_psts_explained) / len(all_psts_known)
statistics.mean(explainability['expl_all_psts'])
statistics.mean(explainability['expl_known_psts'])

all_tprot_explained = set(enz_sub.loc[enz_sub['Pst'].isin(all_psts_explained), 'Tprot'])
all_tprot_pert = set(obs_pert['tprot'])
all_tprot_known = set(all_tprot_pert) & set(enz_sub_univ['tprot'])

len(all_tprot_explained) / len(all_tprot_pert)
len(all_tprot_explained) / len(all_tprot_known)

venn.venn3_unweighted([all_psts_pert, set(enz_sub_univ['pst']), all_psts_explained], set_labels=('psts_pert', 'enz_sub', 'psts_explained'))
venn.venn3_unweighted([all_tprot_pert, set(enz_sub_univ['tprot']), all_tprot_explained], set_labels=('tprot_pert', 'enz_tprot', 'tprot_explained'))
