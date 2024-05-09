from pyScripts.queryProlog import queryProlog
import pandas as pd

# selecting cell line
cline = 'HL60'  # HL60, MCF7, NTERA2

# mainconsultList is a list with files to consult, such as:
mainConsultList = ['rule_dirt.pl',
                   'facts/lint_approach/lint_' + cline + '_mp.pl']

# PROLOG QUERY
df = pd.DataFrame(queryProlog('direct_target(K, Pst, Tprot)', mainConsultList))

# extract enzyme -> target protein relationships (PPIs)
df_tprot = df[['K', 'Tprot']].drop_duplicates()

# export to csv with:
df.to_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/lint_approach/dirts_' + cline + '_mp.csv', index=False)
df_tprot.to_csv('/Users/magdalena/Desktop/bezzlab/logic_things/facts_csv/lint_approach/dirts_tprot_' + cline + '_mp.csv', index=False)

# write Prolog file
with open('facts/lint_approach/dirt_' + cline + '_mp.pl', 'w') as file:
    for index, row in df.iterrows():
        var1 = "dirt('{}', '{}', '{}').".format(row['K'], row['Pst'], row['Tprot'])
        file.write(var1 + '\n')
