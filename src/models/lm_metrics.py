from igraph import Graph
from sklearn.metrics import confusion_matrix, classification_report


def explain_psts(pert_in, edges_df, psts_pert, psts_known):

    # build graph
    g = Graph.DataFrame(edges_df, directed=True, use_vids=False)

    # identify all nodes (psts and tprots) that can potentially be reached from input (perturbagen)
    nodes_reached = g.vs[g.subcomponent(pert_in, mode="out")]["name"]

    # intersect between perturbed and 'reachable' psts
    psts_explained = list(set(psts_pert) & set(nodes_reached)) # returned as list

    # explainability stats
    # compare to ALL perturbed psts
    expl_all_psts = len(psts_explained) / len(psts_pert)
    # compare to perturbed psts with KNOWN enzyme-substrate relationship(s)
    expl_known_psts = len(psts_explained) / len(set(psts_pert) & set(psts_known))
    explainability = [pert_in, len(psts_explained),
                      len(psts_pert), len(set(psts_pert) & set(psts_known)),
                      expl_all_psts, expl_known_psts]

    return psts_explained, explainability


def add_psts_to_k(kins, edges_df, psts_pred):

    k_sink_expl = []

    for k in kins:
        psts_explained, explainability = explain_psts(k, edges_df, psts_pred, psts_pred)
        k_sink_expl.append([k, len(psts_explained), explainability[5]])

    return k_sink_expl


def predict_kins(psts_pred, edges_df, k_sink):

    # build graph
    g = Graph.DataFrame(edges_df, directed=True)

    # empty list of predicted kinases
    k_pred = []

    # identify all kinases that can be reached from at least one perturbed pst
    for k in k_sink:
        asp = g.get_all_shortest_paths(k, psts_pred, mode='out')
        if len(asp) > 0:
            k_pred.append(k)

    return k_pred


def get_pred_metrics(k_pred, pk_true):

    pk_pred = pk_true.fromkeys(pk_true, 0)

    for k in k_pred:
        if k in pk_pred.keys():
            pk_pred[k] = 1
        else:
            pk_pred.update({k: 1})
            pk_true.update({k: 0})

    list(pk_pred.values())
    list(pk_true.values())

    tn, fp, fn, tp = confusion_matrix(list(pk_true.values()), list(pk_pred.values())).ravel()
    cm = {'tn': tn, 'fp': fp, 'fn': fn, 'tp': tp}
    report = classification_report(list(pk_true.values()), list(pk_pred.values()), output_dict=True)

    return cm, report


# NOTE
# v_df = g.get_vertex_dataframe()
# asp = [[g.vs[i]["name"] for i in path] for path in asp]
