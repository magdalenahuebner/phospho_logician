import pandas as pd
import numpy as np
from pyScripts.lm_metrics import get_pred_metrics, add_psts_to_k

# INPUT
obs = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/ctamdb_data/ctamdb_dpoa_MCF7.tsv', sep='\t')
# EDGES
enz_tprot = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/results_lm/enz_tprot_fges_30_bs10_AND_lm_MCF7.csv')
enz_sub = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/results_lm/enz_sub_fges_30_bs10_AND_lm_MCF7.csv')
# UNIVERSE
k_univ = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/r_things/kpas_data/kinases_final.csv')
# VALIDATION
pert_kin = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pert_kin_rel.csv')


# perturbagens that can be validated
pk_true_df = pd.crosstab(index=pert_kin['kinase'], columns = pert_kin['perturbagen'])
perturbagens = pk_true_df.columns.to_list()

# subsetting data
# perturbed (affected) psts IN
obs_pert = obs.loc[(abs(obs['fc']) > 0) &
                   (obs['sid_score'] <= 0.05) &
                   (obs['perturbagen'].isin(perturbagens))].reset_index(drop=True)


# PREDICT
pred_perform = []
pk_pred = dict()
k_predict = []

for pert_in in perturbagens:

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
                             (obs_pert['pst'].isin(edges_df['to'])), 'pst'].tolist()

    # add psts to k_sinks
    k_sink_expl = add_psts_to_k(k_sink, edges_df, psts_pred)

    # kinases actually inhibited by perturbagen (i.e. true values)
    pk_true = dict(pk_true_df[pert_in])

    # determine confusion matrix and accuracy score
    cm, report = get_pred_metrics(k_sink, pk_true)

    # add results to dataframe
    pred_perform.append([pert_in, len(psts_pred), len(k_sink),
                         cm['tp'], cm['fp'], cm['tn'], cm['fn'], report['accuracy'],
                         report['1']['precision'], report['1']['recall'], report['1']['f1-score'], report['1']['support'],
                         report['0']['precision'], report['0']['recall'], report['0']['f1-score'], report['0']['support']])
    pk_pred.update({pert_in: k_sink_expl})
    k_predict = list(np.unique(k_predict + k_sink))

# model performance
pred_perform = pd.DataFrame(pred_perform, columns=['pert_in', 'psts_pred', 'k_pred', 'TP', 'FP', 'TN', 'FN', 'accuracy',
                                                   'precision_1', 'recall_1', 'f1-score_1', 'support_1',
                                                   'precision_0', 'recall_0', 'f1-score_0', 'support_0'])

# predicted pk_rel
pk_pred_df = []
[[pk_pred_df.append([p] + k) for k in pk_pred[p]] for p in pk_pred.keys()]
pk_pred_df = pd.DataFrame(pk_pred_df, columns=['pert_in', 'k_pred', 'psts_explained', 'expl_known_psts'])

# export to csv with:
pred_perform.to_csv('/Users/magdalena/Desktop/bezzlab/logic_things/results_lm/pred_fges_AND_lm_MCF7.csv', float_format='%.3f', index=False)
pk_pred_df.to_csv('/Users/magdalena/Desktop/bezzlab/logic_things/results_lm/pred_pk_fges_AND_lm_MCF7.csv', float_format='%.3f', index=False)


# HEATMAP