import pandas as pd
import seaborn as sns
import scipy
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import math



# function to detect stripeys
#outputs dataframe with a binary stripey call, and a p-value from the Kolmogorov-Smirnov test comparing the timepoint with each of its neighbors. From this information, a stripey will be determined in the binary stripey call, taking into account potential sequential stripeys.
#These p-values will be very small, since each timepoint does have a unique distribution, so a standard 0.05 cutoff will not mean anything. However, the p-values can be compared to determine if a given timepoint is more different from its neighbors.
#examples: In Tina's wildtype RNA-seq data analyzed by RSEM, non-stripey timepoints had a p-value ranging from 0.01 to ~1*10^-200, and each of the stripeys had a p-value of 0. In the same data analyzed by cufflinks, the stripey threshhold that matched visual inspection was 10^-43. For Chun Yi's K562 data, the stripeys both had p-values of 0.


def z_score_normalize(df):
    zscore = scipy.stats.zscore(df, axis=1)
    z_df = pd.DataFrame(zscore, index=df.index, columns=df.columns)
    z_df= z_df.dropna()
    return z_df

def KS_dist_matrix(df):
    z_df = z_score_normalize(df)
    results_matrix = pd.DataFrame(columns = df.columns, index=df.columns)
    for timepoint1 in df.columns:
        for timepoint2 in df.columns:
            results_matrix.loc[timepoint1, timepoint2] = scipy.stats.ks_2samp(z_df[timepoint1].values, z_df[timepoint2].values)[0]
    return results_matrix

def WD_matrix(df):
    z_df = z_score_normalize(df)
    results_matrix = pd.DataFrame(columns = df.columns, index=df.columns)
    for timepoint1 in df.columns:
        for timepoint2 in df.columns:
            results_matrix.loc[timepoint1, timepoint2] = scipy.stats.wasserstein_distance(z_df[timepoint1].values, z_df[timepoint2].values)
    return results_matrix

# function to detect stripeys
def stripey_detector(df, dataset_name = None, comp_stat = True, threshold_stat = 0.3, threshold_pval=0, return_dataframe = False, return_threshplot = True):
        z_df = z_score_normalize(df)
        stripe = pd.DataFrame(index = z_df.columns, columns = ["stripe", "comparison timepoint", "p-value", "statistic"])
        Stripey_overload = False
        for i in range(0, len(z_df.columns)):
            if i ==0:
            ### assume first timepoint is not a Stripey
                stripe.iloc[0]=[0, "none", "none", "none"]
            else:
                if stripe.iloc[i-1][0] == 0:
                ### if previous timepoint is not a stripey, compare to left neighbor
                    stat_ln, pval_ln  = scipy.stats.ks_2samp(z_df.iloc[:,i].values, z_df.iloc[:,i-1].values)
                    if comp_stat == False:
                        test_val = pval_ln
                        threshold = threshold_pval
                    else:
                        test_val = stat_ln *-1
                        threshold = threshold_stat *-1
                        
                    if test_val>threshold:
                    # if above threshold, not a stripey
                        stripe.iloc[i]=[0, z_df.columns[i-1], pval_ln, stat_ln]

                    elif test_val <= threshold:
                    # if below threshold, is a stripey
                        stripe.iloc[i] = [1, z_df.columns[i-1], pval_ln, stat_ln]
                    
                elif stripe.iloc[i-1][0] == 1:
                ### if previous timepont was a stripey
                    done = False
                    for j in range(2, 5):
                    # go back through timepoints up to -5
                        if i-j >=0 and done == False:
                            if stripe.iloc[i-j][0]==0:
                            # if not a stripey compare to that previous neighbor. Set done to True once a non-stripey neighbor is found
                                done = True
                                stat_ls, pval_ls = scipy.stats.ks_2samp(z_df.iloc[:,i].values, z_df.iloc[:,i-j].values)
                                if comp_stat == False:
                                    test_val = pval_ls
                                    threshold = threshold_pval
                                else:
                                    test_val = stat_ls *-1
                                    threshold = threshold_stat *-1
                                
                                if test_val > threshold:
                                # if above threshold, not a stripey
                                    stripe.iloc[i]=[0, z_df.columns[i-j], pval_ls, stat_ls]
                                elif test_val <= threshold:
                                # if below, is a stripey
                                    stripe.iloc[i]=[1, z_df.columns[i-j], pval_ls, stat_ls]
                    if done == False:
                    # if no non-stripey neighbor is found in the 5 left neighbors, re-run it but compare back to front instead of front to back
                        Stripey_overload = True
                
                        
        if Stripey_overload == True:
            print("too many stripeys in forward direction - comparing in reverse direction")
            stripe = pd.DataFrame(index = z_df.columns, columns = ["stripe", "comparison timepoint", "p-value", "statistic"])

            for i in range(0, len(z_df.columns)):
                if i == 0:
                    stripe.iloc[len(z_df.columns)-1]=[0, "none", "none", "none"]
                else:
                    if stripe.iloc[len(z_df.columns)-i][0] == 0:
                        stat_ln, pval_ln  = scipy.stats.ks_2samp(z_df.iloc[:,len(z_df.columns)-1-i].values, z_df.iloc[:,len(z_df.columns)-i].values)
                        if comp_stat == False:
                                test_val = pval_ln
                                threshold = threshold_pval
                        else:
                                test_val = stat_ln *-1
                                threshold = threshold_stat *-1
                        if test_val>threshold:
                            stripe.iloc[len(z_df.columns)-1-i]=[0, z_df.columns[len(z_df.columns)-i], pval_ln, stat_ln]

                        elif test_val <= threshold:
                            stripe.iloc[len(z_df.columns)-1-i] = [1, z_df.columns[len(z_df.columns)-i], pval_ln, stat_ln]
                    
                    elif stripe.iloc[len(z_df.columns)-i][0] == 1:
                        done = False
                        for j in range(2, 5):

                            if i-j >=0 and done == False:
                                if stripe.iloc[len(z_df.columns)-i-1+j][0]==0:
                                    done = True
                                    stat_ls, pval_ls = scipy.stats.ks_2samp(z_df.iloc[:,len(z_df.columns)-1-i].values, z_df.iloc[:,len(z_df.columns)-1-i+j].values)
                                    if comp_stat == False:
                                        test_val = pval_ls
                                        threshold = threshold_pval
                                    else:
                                        test_val = stat_ls *-1
                                        threshold = threshold_stat *-1
                                    if test_val > threshold:
                                        stripe.iloc[len(z_df.columns)-i-1]=[0, z_df.columns[len(z_df.columns)-1-i+j], pval_ls, stat_ls]
                                    elif test_val <= threshold:
                                        stripe.iloc[len(z_df.columns)-i-1]=[1, z_df.columns[len(z_df.columns)-1-i+j], pval_ls, stat_ls]
                        if done == False:
                            Stripey_overload = True
                        
            
   ## return visualization
        if comp_stat == False:
            comp_type = "p-value"
        else:
            comp_type = "statistic"
            if return_threshplot == True:
                if Stripey_overload==True:
                    sns.distplot(stripe.iloc[0:len(stripe)-1]["statistic"].values, bins = len(stripe)-1)
                else:
                    sns.distplot(stripe.iloc[1:len(stripe)]["statistic"].values, bins = len(stripe)-1)

                plt.axvline(threshold_stat, 0, color = 'grey', linestyle='--' )
                plt.xlabel("KS Distance Statistic")
                plt.ylabel("Count")
                if dataset_name == None:
                    title_name = "Distribution of Kologomorv-Smirnov Distance"
                else:
                    title_name = "Distribution of Kologomorv-Smirnov Distance - " + dataset_name
                plt.title(title_name)
                plt.legend(["Threshold"], loc=0)



     # print all stripe time points
        print("Thresholding using", comp_type)
        stripes = []
        for i in range(0, len(stripe)):
            if stripe.iloc[i][0]==1:
                stripes.append(stripe.index[i])
        print("stripes: ", stripes)
        if stripes == []:
            print("No stripeys at this threshhold.")
                                
        if return_dataframe == True:   
            return stripe

## function to visualize distance matrix 
def viz_stripey_dist_matrix(df, dataset_name = None, KS_dist = True):
    if KS_dist == True:
        matrix = KS_dist_matrix(df)
        type_dist = "Kologomorov-Smirnov Distance"
    else:
        matrix = WD_matrix(df)
        type_dist = "Wassserstein Distance"
    matrix = matrix.apply(pd.to_numeric)
    
    sns.heatmap(matrix)
    if dataset_name == None:
        title_name = type_dist + " Matrix"
    else:
        title_name = type_dist + " Matrix - " + dataset_name
    plt.xlabel("Timepoint")
    plt.ylabel("Timepoint")
    plt.title(title_name)

## function to interpolate stripeys
def interpolate_timepoints(df, timepoints, method_option = "pchip"):
    df_int = df.copy()
    for tp in timepoints:
        df_int[tp] = np.nan
    columns_in_df  = list(df.columns)
    df_int.columns = timepoints_numeric = [int(i) for i in columns_in_df]
    df_interpolated = df_int.interpolate(method=method_option, axis =1)
    return df_interpolated

#function to quantial normalize pandas df
#function from Rob Moseley
def qn_normalize(df):
    """
    Perform basic quantile normalization on a pandas dataframe
    """
    qn_df = pd.DataFrame(columns=df.columns)
    temp_df = pd.DataFrame(np.tile(df.apply(np.sort, axis=0).mean(axis=1).values, (len(df.columns),1)).transpose(), columns=df.columns)
    for col in df.columns:
        qn_df[col] = df[col].replace(to_replace=pd.DataFrame(temp_df[col].values, index=df.apply(np.sort, axis=0)[col]).groupby(col).mean().to_dict()[0])

    return qn_df
