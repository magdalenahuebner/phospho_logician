def ksea_ztest(ks_db, data, graphics = False, min_sub = 5):

    """
    test_function does blah blah blah.

    :param ks_db: dataframe with 2 columns - 'kpa' and 'pst' (i.e. Enzyme and Phosphosite)
    :param data: experimental data - dataframe with dimension psts x samples; values: fold changes
    :param graphics: True or False - whether you want to generate graphics
    :param min_sub: minimum number of substrates for a kinase to be considered
    :return: describe what it returns
    """ 
    
    import pandas as pd
    import scipy.stats as st
    import numpy as np
    
    if graphics:
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import axes_size
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
        from mpl_toolkits.axes_grid1.axes_size import AxesY, Fraction
        from mpl_toolkits.axes_grid1.colorbar import colorbar
        import seaborn as sns
    

    # Function used to convert kinase z-scores to corresponding p-values.
    # It is assmued z-scores are normally distributed.
    def get_pvalue(z):
        if z < 0:
            p=st.norm.cdf(z)
            return p
        else:
            dist=st.norm.cdf(z)
            p=1.0 - dist
            return p
        

    # Function to calculate the Z-score for each kinase.
    def get_zscore(mean_kpa, mean_all, sub_num, sd):
        z = (mean_kpa - mean_all) * sub_num**(1/2) / sd
        return z
    
        
    # Function that takes in 2 dataframes: one with kinase/phosphatase-substrate pairs and another with experimental data.
    # It returns a dataframe with z-scores and p-values for each kinase/phosphatase.
    # TODO: implement case where NA values are present in the data
    # TODO: change all_mean and all_std to be calculated per sample
    def run_ksea(enz_sub, data):
        # one hot encode enzsub dataframe to get a matrix which is 1 if pst is a substrate of kpa and 0 otherwise
        enz_sub = enz_sub.assign(val=1).pivot_table(columns='pst',index='kpa',values='val',fill_value=0)
        # calculate the number of substrates for each kinase/phosphatase
        n_sub_data = enz_sub.sum(axis=1)
        # apply matrix multiplication to get the mean fold change for each kinase/phosphatase per sample
        mlog2fc = enz_sub.dot(data).div(n_sub_data, axis=0)
        # calculate the mean fold change and standard deviation for all phosphosites and all samples
        all_mean = data.mean().mean()
        all_std = data.stack().std()
        # calculate the z-score for each kinase/phosphatase and sample
        z_score = mlog2fc.apply(lambda x: get_zscore(x, 0, n_sub_data, 1))
        # calculate the p-value for each kinase/phosphatase and sample
        z_pval = z_score.applymap(lambda x: get_pvalue(x))
        # create a dataframe with kpa name, sample name, n_sub_data, mlog2fc, z_scores and p_values
        ksea_df = pd.DataFrame(z_pval.stack()).reset_index()
        ksea_df.columns = ['kpa', 'pert', 'p_val']
        ksea_df['n_sub_data'] = ksea_df['kpa'].map(n_sub_data)
        ksea_df['mlog2fc'] = mlog2fc.stack().values.tolist()
        ksea_df['z_score'] = z_score.stack().values.tolist()

        return ksea_df
    

    # function that only retains pst if present in both ks_db and data
    def remove_pst(ks_db, data):
        # get pst names from ks_db
        ks_pst = ks_db['pst'].unique()
        # get pst names from data
        data_pst = data.index.tolist()
        # get pst names that are present in both ks_db and data
        pst = [i for i in ks_pst if i in data_pst]
        # only keep pst that are present in both ks_db and data
        ks_db = ks_db[ks_db['pst'].isin(pst)]
        # only keep pst that are present in both ks_db and data
        data = data[data.index.isin(pst)]
        
        return ks_db, data


    # function that calculates n_substrates in ks_db and adds it to ksea_df
    def add_n_sub(ksea_df, ks_db):
        # calculate the number of substrates for each kinase/phosphatase and store in named vector
        n_sub_db = ks_db.groupby('kpa')['pst'].nunique()
        # add n_sub to ksea_df
        ksea_df['n_sub_db'] = ksea_df['kpa'].map(n_sub_db)
        
        return ksea_df
    

    ks_data, data = remove_pst(ks_db, data)
    ksea_df = run_ksea(ks_data, data)
    ksea_df = add_n_sub(ksea_df, ks_db)

    # filter ksea_df to only include kinases/phosphatases where n_sub_data >= min_sub
    ksea_df = ksea_df[ksea_df['n_sub_data'] >= min_sub]
    # make z_score dataframe with kpa as index and sample as columns
    heatmap_df = ksea_df.pivot(index='kpa', columns='pert', values='z_score')
    # make p_value dataframe with kpa as index and sample as columns
    pval_df = ksea_df.pivot(index='kpa', columns='pert', values='p_val')
    # concatenate pval_df into list where each element is a list of p-values for each sample
    pval_list = [pval_df.iloc[:, i].tolist() for i in range(len(pval_df.columns))]
    # p-values for the heatmap annotation are extracted from a nested list into a flat list.
    pvalues=[]
    for entry in pval_list:
        for pval in entry:
            pvalues.append(pval)
            
    
    # Heatmap only generated if the user chose to produce graphics during file upload.
    # TODO: improve colour palette
    if graphics:
        # Set the margins and square height for a single category.
        topmargin = 0.1 #inches
        bottommargin = 0.1 #inches
        categorysize = 0.35 # inches
        # Number of kinases identified.
        n=heatmap_df.shape[0]
        leftmargin = 0.1
        rightmargin = 0.1
        catsize = 0.5
        # Number of conditions (e.g. cell lines).
        m=heatmap_df.shape[1]

        # Parameters for color bar.
        aspect = n
        pad_fraction = 0.7

        # Calculate a dynamic figure height.
        figheight = topmargin + bottommargin + (n+1)*categorysize

        # Calculate a dynamic figure width.
        figwidth = leftmargin + rightmargin + (m+1)*catsize

        fig, ax = plt.subplots(figsize=(figwidth, figheight))

        # Format the axes.
        ax.xaxis.set_ticks_position('top')
        plt.yticks(fontsize=6)
        plt.xticks(fontsize=6)

        # Plot the heatmap.
        ax = sns.heatmap(heatmap_df, cmap='vlag', annot=True, fmt=".1f", annot_kws={'size':5}, cbar=False, linewidths=0.3, linecolor='white')

        # Format the colour bar dynamically.
        ax_div = make_axes_locatable(ax)
        width = axes_size.AxesY(ax, aspect=1./aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax = ax_div.append_axes('right', size = width, pad = pad)
        cb=plt.colorbar(ax.get_children()[0], cax = cax, orientation = 'vertical')
        cax.yaxis.set_ticks_position('right')
        cb.ax.tick_params(labelsize=6)
        cb.set_label('Z-Score', fontsize=6, labelpad=7)
        cb.outline.set_visible(False)

        #Remove y-axis label.
        ax.yaxis.set_label_text("")

        # Rotate the axis labels.
        for item in ax.get_yticklabels():
            item.set_rotation(0)

        for item in ax.get_xticklabels():
            item.set_rotation(90)

        # Annotate statistically significant scores with asterisks.
        # * for p < 0.05 and ** for p < 0.01.
        counter=-1
        for text in ax.texts:
            counter+=1
            if pvalues[counter] < 0.05 and pvalues[counter] >= 0.01:
                text.set_weight('bold')
                text.set_text(text.get_text() + "*")
            elif pvalues[counter] < 0.01:
                text.set_weight('bold')
                text.set_text(text.get_text() + "**")

        # Add a title.
        ax.set_title('Kinase/Phosphatase Activity', fontsize=8, pad=10)
        
        # save fig
        plt.savefig('ksea_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()

    return ksea_df