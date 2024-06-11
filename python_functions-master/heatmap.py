import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import matplotlib



#make colormap haase lab color scheme
norm = matplotlib.colors.Normalize(-1.5,1.5)
colors = [[norm(-1.5), "cyan"],
          [norm(0), "black"],
         [norm(1.5), "yellow"]]
haase = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)




def get_closest_column_from_period(dataset_df, period):
    '''
    Returns the column from the supplied dataset that is closest to the supplied period.

    Parameters
    ----------
    dataset_df : pandas.DataFrame
         time series gene expression dataset, where rows are genes and columns are time points
    period : integer
        period which to find the closest timepoint in dataset_df to

    Returns
    -------
    closest_column_int : integer
        the timepoint in dataset_df's columns which the period is closest to

    '''
    columns_in_df = list(dataset_df.columns)
    timepoints_numeric = [int(i) for i in columns_in_df]

    closest_column_int = timepoints_numeric[min(range(len(timepoints_numeric)), key = lambda i: abs(timepoints_numeric[i]-period))]
    return closest_column_int



## function to plot data in a heat map, sorted based off of the maximum in the first period, and normalized by z-score
def heatmap_max(data, gene_list, first_period, z_score_norm = True, first_period_name = True, yticks = False, cbar_bool = True, axis = None, cmap = haase, scale_min = -1.5, scale_max = 1.5):
    # data = data frame with gene names as the index
    # gene list = list of gene names to plot
    # first period = integer value of the last column of the first period

    # yticks: default False doesn't plot ytick labels
    # cbar_bool: default True does include color bar
    # axis: default None doesn't specify an axis
    data = data.loc[gene_list]
    #use .reindex?
    if z_score_norm == True:
        norm_data = scipy.stats.zscore(data, axis=1)
    else:
        norm_data  = data
    norm_data_df = pd.DataFrame(norm_data, index=data.index, columns=data.columns)

    if first_period_name:
        first_period = get_closest_column_from_period(data, first_period)
        first_period_index = data.columns.get_loc(str(first_period))
    else: 
        first_period_index = first_period


    norm_data_1stper = norm_data_df.iloc[:, 0:first_period_index]
    max_time = norm_data_1stper.idxmax(axis=1)
    norm_data_df["max"] = max_time
    norm_data_df["max"] = pd.to_numeric(norm_data_df["max"])
    norm_data_df = norm_data_df.sort_values(by="max", axis=0)
    norm_data_df = norm_data_df.drop(columns=['max'])
    order = norm_data_df.index
    s = sns.heatmap(norm_data_df, cmap=cmap, vmin=scale_min, vmax=scale_max, yticklabels = yticks, cbar = cbar_bool, ax = axis)
    return order


## function to plot data in a heat map, sorted based off of the supplied order, and normalized by z-score
def heatmap_order(data, order, z_score_norm = True, yticks = False, cbar_bool = True, axis = None, scale_min = -1.5, scale_max = 1.5):
    # data = data frame with gene names as the index
    #order = list of genes to plot in desired order

    # yticks: default False doesn't plot ytick labels
    # cbar_bool: default True does include color bar
    # axis: default None doesn't specify an axis
    data = data.loc[order]
    if z_score_norm == True:
        norm_data = scipy.stats.zscore(data, axis=1)
    else:
        norm_data = data
    norm_data_df = pd.DataFrame(norm_data, index=data.index, columns=data.columns)
    order2 = norm_data_df.index
    s = sns.heatmap(norm_data_df, cmap=haase, vmin=scale_min, vmax=scale_max, yticklabels = yticks, cbar= cbar_bool,ax = axis)
    return order2

## function to plot data in a heat map, sorted based off of the left edge, and normalized by z-score
def heatmap_LE(data,gene_list, first_period,z_score_norm = True, first_period_name =True,  yticks = False, cbar_bool = True, axis = None,  scale_min = -1.5, scale_max = 1.5):
    #adapted from R code orderHeatmapsLE (from R function autoorder)

    # data = data frame with gene names as the index
    # gene list = list of gene names to plot
    # first period = integer value of the last column of the first period

    # yticks: default False doesn't plot ytick labels
    # cbar_bool: default True does include color bar
    # axis: default None doesn't specify an axis
    data = data.loc[gene_list]
    if z_score_norm == True:
        norm_data = scipy.stats.zscore(data, axis=1)
    else:
        norm_data = data
    norm_data_df = pd.DataFrame(norm_data, index=data.index, columns=data.columns)

    if first_period_name:
        first_period = get_closest_column_from_period(data, first_period)
        first_period_index = data.columns.get_loc(str(first_period))
    else: 
        first_period_index = first_period

    norm_data_1stperiod = norm_data_df.iloc[:, 0:first_period_index]    
    LE = pd.DataFrame(index=norm_data_1stperiod.index)
    LE["order"] = 0
    for i in range(0, len(norm_data_1stperiod.index)):
        temp = norm_data_1stperiod.iloc[i]
        sPos = -1
        maxLength = 0
        for j in range(0, len(temp)):
            if temp.iloc[j] > 1:
                if sPos == -1:
                    sPos = j
                    curLength = 0
                curLength = curLength + 1
            if temp.iloc[j] < 1:
                if sPos != -1:
                    if curLength > maxLength:
                        maxLength = curLength
                        LE.iloc[i] = sPos
        if sPos != -1:
            if curLength > maxLength:
                maxLength = curLength
                LE.iloc[i] = sPos
    norm_data_df["order"] = LE["order"]
    norm_data_df = norm_data_df.sort_values(by="order", axis=0)
    norm_data_df = norm_data_df.drop(columns=['order'])
    order = norm_data_df.index
    s = sns.heatmap(norm_data_df, cmap=haase, vmin=scale_min, vmax=scale_max, yticklabels = yticks, cbar = cbar_bool, ax = axis )

    return order



# def heatmap_LE_max(data, gene_list, first_period, yticks = False, cbar_bool = True):
#     # data = data frame with gene names as the index
#     # gene list = list of gene names to plot
#     # first period = integer value of the last column of the first period
#     data = data.loc[gene_list]
#     z_pyjtk = scipy.stats.zscore(data, axis=1)
#     z_pyjtk_df = pd.DataFrame(z_pyjtk, index=data.index, columns=data.columns)
#     z_pyjtk_1stperiod = z_pyjtk_df.iloc[:, 0:first_period]
#     max_time = z_pyjtk_1stperiod.idxmax(axis=1)
#     z_pyjtk_df["max"] = max_time
#     # z_pyjtk_df["max"] = pd.to_numeric(z_pyjtk_df["max"])
#     z_pyjtk_df = z_pyjtk_df.sort_values(by="max", axis=0)
#     z_pyjtk_df = z_pyjtk_df.drop(columns=['max'])
#     z_pyjtk_1stperiod = z_pyjtk_df.iloc[:, 0:first_period]
#
#
#     LE = pd.DataFrame(index=z_pyjtk_1stperiod.index)
#     LE["order"] = 0
#     for i in range(0, len(z_pyjtk_1stperiod.index)):
#         temp = z_pyjtk_1stperiod.iloc[i]
#         sPos = -1
#         maxLength = 0
#         for j in range(0, len(temp)):
#             if temp.iloc[j] > 1:
#                 if sPos == -1:
#                     sPos = j
#                     curLength = 0
#                 curLength = curLength + 1
#             if temp.iloc[j] < 1:
#                 if sPos != -1:
#                     if curLength > maxLength:
#                         maxLength = curLength
#                         LE.iloc[i] = sPos
#         if sPos != -1:
#             if curLength > maxLength:
#                 maxLength = curLength
#                 LE.iloc[i] = sPos
#     z_pyjtk_df["order"] = LE["order"]
#     z_pyjtk_df = z_pyjtk_df.sort_values(by="order", axis=0)
#     z_pyjtk_df = z_pyjtk_df.drop(columns=['order'])
#     order = z_pyjtk_df.index
#
#     max_time = z_pyjtk_1stperiod.idxmax(axis=1)
#     z_pyjtk_df["max"] = max_time
#     # z_pyjtk_df["max"] = pd.to_numeric(z_pyjtk_df["max"])
#     z_pyjtk_df = z_pyjtk_df.sort_values(by="max", axis=0)
#     z_pyjtk_df = z_pyjtk_df.drop(columns=['max'])
#     z_pyjtk_1stperiod = z_pyjtk_df.iloc[:, 0:first_period]
#     s = sns.heatmap(z_pyjtk_df, cmap=haase, vmin=-1.5, vmax=1.5, yticklabels = yticks, cbar = cbar_bool)
#
#     return order
