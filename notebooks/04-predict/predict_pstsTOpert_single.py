import pandas as pd
import numpy as np
from pyScripts.lm_metrics import get_pred_metrics

# INPUT
obs = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/observations_HL60.csv')
# EDGES
enz_tprot = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/enz_tprot_strong_HL60.csv')
enz_sub = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/enz_sub_strong_HL60.csv')
# UNIVERSE
k_univ = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/r_things/kpas_data/kinases_final.csv')
# VALIDATION
pert_kin = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pert_kin_rel.csv')


# perturbagens that can be validated
pk_true_df = pd.crosstab(index=pert_kin['kinase'], columns = pert_kin['perturbagen'])
perturbagens = pk_true_df.columns.to_list()
pert_in = 'Torin'  # select one

# subsetting data
# perturbed (affected) psts IN
obs_pert = obs.loc[(abs(obs['fold_change']) > 0) &
                   (obs['p_value'] <= 0.05) &
                   (obs['perturbagen'] == pert_in)].reset_index(drop=True)


# PREDICT
# edges_df
enz_tprot_edges = np.array(enz_tprot.loc[enz_tprot['Pert'] == pert_in, ['Kpa', 'Tprot']])
enz_sub_edges = np.array(enz_sub.loc[enz_sub['Pert'] == pert_in, ['Kpa', 'Pst']])
edges_df = pd.DataFrame(np.concatenate((enz_tprot_edges, enz_sub_edges)),
                        columns=['from', 'to'])

# sink kinases (i.e. predicted kinases)
nodes_sink = set(edges_df['from']) - set(edges_df['to'])
k_sink = list(set(nodes_sink) & set(k_univ['kinase']))

# psts affected by perturbagen and contained in edges_df (i.e. predictive psts)
psts_pred = obs_pert.loc[(obs_pert['perturbagen'] == pert_in) &
                         (obs_pert['substrate'].isin(edges_df['to'])), 'substrate'].tolist()

# kinases actually inhibited by perturbagen (i.e. true values)
pk_true = dict(pk_true_df[pert_in])

# determine confusion matrix and accuracy score
cm, accuracy = get_pred_metrics(k_sink, pk_true)

# model performance
pred_perform = [pert_in, len(psts_pred), len(k_sink), cm['tp'], cm['fp'], cm['tn'], cm['fn'], accuracy]


# HEATMAP
