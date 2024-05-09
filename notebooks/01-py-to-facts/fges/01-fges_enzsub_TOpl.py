import pandas as pd

enz_sub = pd.read_csv('/Users/magdalena/Desktop/bezzlab/logic_things/fges_data/ctamdb_data/results/fges_enzsub20_bs10_MCF7.txt', sep='\t')

# enz_trot
enz_sub['tprot'] = enz_sub['pst'].str.replace(r'\(.*$', "")

# write Prolog file
with open('facts/fges/fges_enzsub20_bs10_MCF7.pl', 'w') as file:
    for index, row in enz_sub.iterrows():
        var1 = "fges_enzsub('{}', '{}', '{}').".format(row['kpa'], row['pst'], row['tprot'])
        file.write(var1 + '\n')
