# Function that splits a list in n elements
def split_it(seq, n):

    avg = len(seq) / float(n)
    lst = []
    count = 0.0

    while count < len(seq):
        lst.append(seq[int(last):int(last + avg)])
        count += avg

    return lst

################################################################################

# Function that impute missing values (NaN) with the minimum value of the row
def min_global_imputation(df):

    # Iterate through rows
    for index, row in df.iterrows():
        # Get the minimum value of the row
        min_value = df.loc[index, :].min()
        # Replace NaNs with the minimum value
        df.loc[index, :].replace(np.nan, min_value, inplace = True)

    return df

################################################################################

# Function that imput missing values (NaN) with the minimum values of a
# treatment arm (while working with replicates)
def min_local_imputation(df):

    # Splits the columns in the respective treatment arms
    cols = split_it(list(df.columns), 8)

    # Iterate through rows
    for index, row in df.iterrows():
        # Iterate through the treatments
        for treatment in cols:
            count = 0
            # Iterate through the replicates of a treatment arm
            for sample in treatment:
                # Evaluates if the value is a NaN and score it in the count
                if np.isnan(df.loc[index, sample]):
                    count += 1

            # If there is one or two missing values in a treatment the NaN will be replaced with the minimum value of the arm
            if count == 1 or count == 2:
                min_value = df.loc[index, treatment].min()
                for samp in treatment:
                    if np.isnan(df.loc[index, samp]):
                        df.loc[index, samp] = min_value
                        continue

                continue

    return df

################################################################################

# Calculate multiple testing anova test and corrects their p-values
def anova_test(df, type, correction):

    # Creates empty df with the columns that are the output of Pingouin anova functions
    anova = pd.DataFrame(columns = ['Protein', 'ddof1', 'ddof2', 'F', 'p-unc', 'np2'])

    # Iterate through unique proteins
    for protein in df_stats['Accession'].unique():

        # Create temporal df for each protein
        df_fil = df[df['Accession'] == protein]

        # Compute type of anova for each protein
        if type == 'anova':
            aov = pg.anova(dv = 'Intensity', between = 'Experiment', data = df_fil)
        elif type == 'welch':
            aov = pg.welch_anova(dv = 'Intensity', between = 'Experiment', data = df_fil)
        else:
            print('Anova type error')


        # Append the anova results to the new df
        anova = anova.append({'Protein': protein,
                             'ddof1': aov['ddof1'][0],
                             'ddof2': aov['ddof2'][0],
                             'F': aov['F'][0],
                             'p-unc': aov['p-unc'][0],
                             'np2': aov['np2'][0]}, ignore_index = True)

    # Compute the adjustes p-values and append them to the new df
    pvals = anova.loc[:, 'p-unc'].values
    anova['p-adj'] = pg.multicomp(pvals, method = correction)[1]

    return anova

################################################################################
