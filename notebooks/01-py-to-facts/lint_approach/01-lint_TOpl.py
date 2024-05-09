import pandas as pd

# selecting cell line
cline = 'MCF7'  # HL60, MCF7, NTERA2

# load known perturbagen-kinase relationships
pk = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/pert_kin_rel.csv')
# load observational data (has to be repeated for every cell line individually)
obs = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/ctamdb_data/ctamdb_dpoa_' + cline + '.tsv', sep='\t')
# load negatome
ngt = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/r_things/negatome/negatome.csv')

pks = []  # perturbagen, kinase, substrate
pks_neg = []  # negative statements

for i, row in pk.iterrows():
    print(i)
    p = row['perturbagen']
    k = row['kinase']
    subs = obs.loc[obs['perturbagen'] == p]
    for j, s in subs.iterrows():
        if s['fc'] != 0 and s['sid_score'] <= 0.05:
            pks.append([p, k, s['pst']])
        else:
            pks_neg.append([p, k, s['pst']])

pks = pd.DataFrame(pks, columns=['perturbagen', 'kinase', 'pst'])
pks_neg = pd.DataFrame(pks_neg, columns=['perturbagen', 'kinase', 'pst'])

lint_pos = pks[['kinase', 'pst']].drop_duplicates()
lint_neg = pks_neg[['kinase', 'pst']].drop_duplicates()

# use negative statements to update positive statements
lint_corr = pd.concat([lint_pos, lint_neg, lint_neg]).drop_duplicates(keep=False)

# add also the target protein (the protein that a phosphosite is on) to df
lint_pos['tprot'] = lint_pos['pst'].str.replace(r'\(.*$', "")
lint_corr['tprot'] = lint_corr['pst'].str.replace(r'\(.*$', "")

# use negatome to update statements
ngt.columns = ['kinase', 'tprot']
lint_corr = pd.merge(lint_corr, ngt, how='outer', indicator=True)
lint_corr = lint_corr.loc[lint_corr._merge == 'left_only', ['kinase', 'pst', 'tprot']]

# export df to csv
lint_corr.to_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/lint_approach/lints_' + cline + '_mp.csv', index=False)
lint_pos.to_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/lint_approach/lints_ps_' + cline + '_mp.csv', index=False)

# write Prolog file
with open('facts/lint_approach/lint_' + cline + '_mp.pl', 'w') as file:
    for index, row in lint_corr.iterrows():
        var1 = "us_of('{}', '{}', '{}').".format(row['kinase'], row['pst'], row['tprot'])
        file.write(var1 + '\n')
