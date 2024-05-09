import pandas as pd
import numpy as np
from pyScripts.lm_metrics import explain_psts

# INPUT
pert_kin = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pert_kin_rel.csv')
# EDGES
enz_sub = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/fges_data/ctamdb_data/results/fges_enzsub20_bs10_MCF7.txt', sep='\t')
# OUTPUT
obs = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/ctamdb_data/ctamdb_dpoa_MCF7.tsv', sep='\t')

# enz_trot
enz_tprot = enz_sub.copy()
enz_tprot['tprot'] = enz_tprot['pst'].str.replace(r'\(.*$', "")
enz_tprot = enz_tprot[['kpa', 'tprot']].drop_duplicates()

# perturbagens IN
perturbagens = pert_kin['perturbagen'].drop_duplicates().tolist()

# subsetting data
# perturbed (affected) psts OUT
obs_pert = obs.loc[(abs(obs['fc']) > 0) &
                   (obs['sid_score'] <= 0.05) &
                   (obs['perturbagen'].isin(perturbagens))].reset_index(drop=True)
# psts KNOWN in enzyme-substrate relationship(s) universe
psts_known = list(set(enz_sub['pst']))


# EXPLAIN
explainability = []
psts_explained = []

# edges_df
pert_kin_edges = np.array(pert_kin[['perturbagen', 'kinase']])
enz_tprot_edges = np.array(enz_tprot)
enz_sub_edges = np.array(enz_sub[['kpa', 'pst']])
edges_df = pd.DataFrame(np.concatenate((pert_kin_edges, enz_tprot_edges, enz_sub_edges)),
                        columns=['from', 'to'])

for pert_in in perturbagens:

    # psts affected by perturbagen
    psts_pert = obs_pert.loc[obs_pert['perturbagen'] == pert_in, 'pst'].tolist()

    # determine explainable phosphosites and explainability
    psts_expl, expl = explain_psts(pert_in, edges_df, psts_pert, psts_known)

    # add results to dataframe
    explainability.append(expl)
    psts_explained = list(np.unique(psts_explained + psts_expl))

explainability = pd.DataFrame(explainability, columns=['pert_in', 'psts_explained', 'psts_pert', 'psts_known', 'expl_all_psts', 'expl_known_psts'])


# export to csv with:
explainability.to_csv('/Users/magdalena/Desktop/bezzlab/logic_things/fges_data/ctamdb_data/results/expl_fges_MCF7.csv', index=False)
