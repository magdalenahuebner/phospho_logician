import pandas as pd
import numpy as np
from pyScripts.lm_metrics import predict_kins, get_pred_metrics

# INPUT
obs = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/observations_HL60.csv')
# EDGES
enz_tprot = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdts_tprot_HL60.csv')
enz_sub = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pdts_HL60.csv')
# UNIVERSE
k_univ = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/r_things/kpas_data/kinases_final.csv')
# VALIDATION
pert_kin = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pert_kin_rel.csv')


# perturbagens that can be validated
pk_true_df = pd.crosstab(index=pert_kin['kinase'], columns = pert_kin['perturbagen'])
perturbagens = pk_true_df.columns.to_list()
k_valid = pk_true_df.index.to_list()

# subsetting data
# perturbed (affected) psts IN
obs_pert = obs.loc[(abs(obs['fold_change']) > 0) &
                   (obs['p_value'] <= 0.05) &
                   (obs['perturbagen'].isin(perturbagens))].reset_index(drop=True)


# PREDICT
pred_perform = []
pk_pred = dict()
k_predict = []

# edges_df
enz_tprot_edges = np.array(enz_tprot)
enz_sub_edges = np.array(enz_sub[['kinase', 'substrate']])
edges_df = pd.DataFrame(np.concatenate((enz_tprot_edges, enz_sub_edges)),
                        columns=['from', 'to'])

# sink kinases
nodes_sink = set(edges_df['from']) - set(edges_df['to'])
k_sink = list(set(nodes_sink) & set(k_univ['kinase']))

for pert_in in perturbagens:

    # psts affected by perturbagen and contained in edges_df (i.e. predictive psts)
    psts_pred = obs_pert.loc[(obs_pert['perturbagen'] == pert_in) &
                             (obs_pert['substrate'].isin(edges_df['to'])), 'substrate'].tolist()

    # kinases actually inhibited by perturbagen (i.e. true values)
    pk_true = dict(pk_true_df[pert_in])

    # determine predicted kinases from sink kinases
    k_pred = predict_kins(psts_pred, edges_df, k_sink)

    # determine confusion matrix and accuracy score
    cm, report = get_pred_metrics(k_pred, pk_true)

    # add results to dataframe
    pred_perform.append([pert_in, len(psts_pred), len(k_pred),
                         cm['tp'], cm['fp'], cm['tn'], cm['fn'], report['accuracy'],
                         report['1']['precision'], report['1']['recall'], report['1']['f1-score'], report['1']['support'],
                         report['0']['precision'], report['0']['recall'], report['0']['f1-score'], report['0']['support']])
    pk_pred.update({pert_in: k_pred})
    k_predict = list(np.unique(k_predict + k_pred))

# model performance
pred_perform = pd.DataFrame(pred_perform, columns=['pert_in', 'psts_pred', 'k_pred', 'TP', 'FP', 'TN', 'FN', 'accuracy',
                                                   'precision_1', 'recall_1', 'f1-score_1', 'support_1',
                                                   'precision_0', 'recall_0', 'f1-score_0', 'support_0'])


# export to csv with:
pred_perform.to_csv('/Users/magdalena/Desktop/bezzlab/logic_things/results_lm/pred_pdts_HL60.csv', float_format='%.3f', index=False)


# HEATMAP (DO IN R)
