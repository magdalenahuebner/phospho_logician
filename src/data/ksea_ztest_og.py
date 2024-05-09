def kseaZtest(kt_db, df):
    global tnum_dic
    import pandas as pd
    import scipy.stats as st
    import numpy as np

    # Function to calculate the Z-score for each kinase or phosphatase.
    def zScore(mean_kpa, mean_all, t_num, sd):
        z = (mean_kpa - mean_all) * t_num ** (1 / 2) / sd
        return z

    # Function used to convert kinase/phosphatase z-scores to corresponding p-values.
    def getpValue(z):
        if z < 0:
            p = st.norm.cdf(z)
            return p
        else:
            dist = st.norm.cdf(z)
            p = 1.0 - dist
            return p

    # User data is parsed as appropriate. Sample conditions (e.g. perturbagens) are extracted from header
    header = df.columns.values.tolist()
    col_length = len(header)
    array = df.to_numpy()

    kpa_dic = {}
    dic = {}
    kt_links = []
    zscore_info = []
    z_columns = ["kpa", "tc_total"]

    ### COLUMN LOOP START ###

    # Columns 1 and onwards represent samples (e.g. perturbagens).
    # For the current column (sample) a set of operations is performed.
    for col in range(1, col_length):
        # Column names for relevant dataframes are created here dynamically.
        # curr_col is current sample/column name.
        curr_col = header[col]
        z_columns.append("tc." + curr_col)
        z_columns.append("mlog2fc." + curr_col)
        z_columns.append("z_score." + curr_col)
        z_columns.append("p_val." + curr_col)

        # Final dictionary contains unique phosphosites and individual log2(FC) values. Some phosphosites were only
        # measured in some of the samples. Since the input dataframe requires a value for each phosphosite/sample
        # pair, some phosphosite/sample pairs might produce nans.
        for n in range(0, len(array)):
            site = array[n][0]
            fc = array[n][col]
            dic[site] = fc

        # For each sample, the mean log2(FC) and standard deviation of all phosphosites in the dataset are calculated
        # here. These values will be used to obtain a z-score for each identified kpa later on.
        all_log2 = []
        for key in dic:
            # Some phosphosites don't have log2(FC) since they weren't measured in a particular sample. They are
            # excluded from calculations.
            if dic[key] == 'nan':
                continue
            else:
                all_log2.append(dic[key])
        all_mean = sum(all_log2) / float(len(all_log2))
        all_std = np.std(all_log2)

        # Each phosphosite in the dictionary is scanned against the known_targets db.
        # If a match is found, relevant information for that phosphosite is retained.
        # Scanning is only done for the first column.
        if col == 1:
            for x in dic:
                print(x)
                kpas = list(kt_db.loc[kt_db['pst'] == x, 'kpa'])
                if kpas == []:
                    continue
                else:
                    for y in kpas:
                        # kt_links will be used to assign the current sample's log2(FCs) to each kpa later on.
                        kt_links.append([y, x, dic[x]])
                        # Following dictonary contains the number of substrate psts for each kinase that were
                        # measured in at least one of the conditions (maximum nr. of potentially measured substrates)
                        tnum_dic = dict(pd.DataFrame(kt_links).groupby([0]).size())
        # Once the first column is passed, new log2(FCs) are removed and/or appended to the original arrays for each sample.
        elif col > 1:
            for s in kt_links:
                s.remove(s[-1])
                s.append(dic[s[1]])

        # A dictionary containing unique kpas and pst log2(FCs) is created.
        # Reset the values in the dictionary for a new column.
        for kpa in kpa_dic:
            kpa_dic[kpa] = []
        # If the same kpa was identified for multiple psts, multiple log2(FCs) are appended to the dictionary values.
        for match in kt_links:
            kpa = match[0]
            log2fc = match[-1]
            if kpa not in kpa_dic:
                kpa_dic[kpa] = [log2fc]
            else:
                kpa_dic[kpa].append(log2fc)

        # The dictionary is used to calculate the number of substrates identified for each kpa.
        # It also calculates the mean log2(FC) across each kpa's targets.
        # The algorithm then computes the z-score.
        # A new array stores kpa, no. of substrates, mean log2(FC), z-score and p-value.
        index = -1
        for key in kpa_dic:
            index += 1
            tnum_total = tnum_dic[key]
            # nans are included in the mapping so that the dict length for every sample is the same (all kpas are
            # included, even if no substrate psts were measured). However, they have to be excluded for calculating z-scores.
            fcs = [elem for elem in kpa_dic[key] if elem != 'nan']
            target_num = len(fcs)
            if target_num == 0:
                if col == 1:
                    zscore_info.append([key, tnum_total, target_num, 'nan', 'nan', 'nan'])
                else:
                    zscore_info[index].append(target_num)
                    zscore_info[index].append('nan')
                    zscore_info[index].append('nan')
                    zscore_info[index].append('nan')
            else:
                kpa_fc_mean = sum(fcs) / float(target_num)
                z_score = zScore(kpa_fc_mean, all_mean, target_num, all_std)
                p_val = getpValue(z_score)
                if col == 1:
                    # KSEA results stored here.
                    zscore_info.append([key, tnum_total, target_num, kpa_fc_mean, z_score, p_val])
                # If the program has gone past the first condition column, kin_fc_mean, z_score and zpval are
                # appended in a repeating manner to the original array.
                else:
                    zscore_info[index].append(target_num)
                    zscore_info[index].append(kpa_fc_mean)
                    zscore_info[index].append(z_score)
                    zscore_info[index].append(p_val)

    ### COLUMN LOOP END ###

    # Score and kinase-substrate relationships DFs are generated.
    zscore_df = pd.DataFrame(zscore_info, columns=z_columns)
    zscore_df = zscore_df.sort_values(by=['kpa'], ignore_index=True)

    return zscore_df


# Function to reformat output dataframe from kseaZtest function
def stackSamples(df):
    import pandas as pd
    header = df.iloc[:, 2:].columns.values.tolist()
    tuples = [tuple(i.split('.', 1)) for i in header]
    midx = pd.MultiIndex.from_tuples(tuples, names=('parameter', 'pert'))
    df_idx = pd.DataFrame(df.iloc[:,2:].to_numpy(), columns=midx, index=df['kpa'])
    df_stack = df_idx.stack()
    df_stack.reset_index(inplace=True, level=1)
    df_stack = pd.merge(df.iloc[:, :2], df_stack, left_on='kpa', right_index=True)
    df_stack = df_stack.reset_index(drop=True)
    
    return df_stack