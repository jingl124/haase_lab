import pandas as pd
import os
from scipy.optimize import curve_fit
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import numpy as np
import structs
import importlib
import graphviz
import matplotlib
import math
import datetime

importlib.reload(structs)

# global dictionary of GeneNode names (strings) and the actual GeneNodes that already exist
gene_nodes = {}

# global dictionary of expression data
exp_dict = {}

# extracting and wrangling data
def extract_expression_data(df):
    '''
    This function extracts the necessary information from the DataFrame and stores it into global dictionary exp_dict.
    exp_dict is formatted like this:
        {
            TF1: {
                target_gene1: [(time, log2_cleaned_ratio), (time, log2_cleaned_ratio), ...],
                target_gene2: [(time, log2_cleaned_ratio), (time, log2_cleaned_ratio), ...],
                ...
            }
            TF2: {
                target_gene1: [(time, log2_cleaned_ratio), (time, log2_cleaned_ratio), ...],
                target_gene2: [(time, log2_cleaned_ratio), (time, log2_cleaned_ratio), ...],
                ...
            }
        }
    Parameters:
    df (pandas DataFrame): contains all data derived from gene expression data file
    '''
    for _, row in df.iterrows():
        tf = row['TF']
        gene = row['GeneName']
        if tf not in exp_dict:
            exp_dict[tf] = {}  # Initialize nested dictionary for tf if not present
        if gene not in exp_dict[tf]:
            exp_dict[tf][gene] = []  # Initialize list for gene if not present
        exp_dict[tf][gene].append([row['time'],
                                   row['log2_cleaned_ratio']])

def sort_genes(path, sep):
    '''
    Returns the list of TFs and targets to be included in GRN based on information found in an input file. 
    Parameters:
    path (string): path of the input file
    sep (string): separator to parse input file

    Returns:
    tfs (list of strings): list of TFs
    targets (list of strings): list of target genes
    '''
    df = pd.read_csv(path, sep=sep)
    tfs = []
    targets = []
    for _, row in df.iterrows():
        gene = row['gene']
        gene_type = row['type']
        if 'target' in gene_type:
            targets.append(gene)
        if 'tf' in gene_type:
            tfs.append(gene)
    return tfs, targets

def filter_df(path, gene_path):
    '''
    Read input files to store data of desired TFs and targets into a pandas DataFrame.

    Parameters: 
    path (string): input file containing all gene expression data of all TFs and targets
    gene_path (string): input file containing TFs and targets

    Returns: 
    df (pandas DataFrame): contains gene expression data of TFs and targets
    '''
    # reading data
    df = pd.read_csv(path, sep='\t')

    # restricting nodes
    tfs, targets = sort_genes(gene_path, ',')
    df = df[df['TF'].isin(tfs) & df['GeneName'].isin(targets)]
    if df.empty:
        raise Exception("DataFrame is empty. Please check the input TFs and target genes.")
    
    # filter columns
    df = df[['TF', 'GeneName', 'time', 'log2_cleaned_ratio']]
    return df

def read_peak_times():
    '''
    Read peak expression times of genes.

    Returns:
    peak_times (pandas DataFrame): peak times of present genes
    '''
    df = pd.read_csv('peak_times.csv')
    peak_times = df[df['Genes'].isin(gene_nodes)]
    return peak_times

# using sigmoidal fit to retrieve gene expression info
def get_sig_info(tf, gene, amp_thresh):
    '''
    Retrieve whether the interaction is activation/inhibition and when it occurs based on sigmoidal curve fit. 
    This is based on the t1/2 of the sigmoidal curve.
    
    Parameters:
    tf (string): name of TF
    gene (string): name of target gene
    amp_thresh (int or float): the minimum amplitude of the sigmoidal curve to be considered activation/inhibition

    Returns:
    time of activation/inhibition (float)
    whether its activation or not (boolean)
    '''
    params = sigmoid_curve_fit(tf, gene) # [L_opt, x0_opt, k_opt, b_opt] or [L1, x01, k1, b1, L2, x02, k2, b2]
    if params is None or len(params) == 0 or abs(params[0]) < amp_thresh:
        return None, None
    return params[1], (params[0] > 0)

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return y

def double_sigmoid(x, L1, x01, k1, b1, L2, x02, k2, b2):
    y = (L1 / (1 + np.exp(-k1 * (x - x01))) + b1) + (L2 / (1 + np.exp(-k2 * (x - x02))) + b2)
    return y

def scale_data(data):
    '''
    Scale the data to a range suitable for fitting.

    Parameters: 
    data (list of floats): expression level data for one TF-target pair

    Returns:
    (list of floats): rescaled version of expression level data
    '''
    data_min = min(data)
    data_max = max(data)
    if data_max == data_min:
        return [0 for _ in data]
    return [(d - data_min) / (data_max - data_min) for d in data]

def sigmoid_curve_fit(tf, gene):
    '''
    Find the sigmoidal curve of best fit (single or double).

    Parameters:
    tf (string): name of TF
    gene (string): name of target gene

    Returns: 
    params (list): sigmoidal curve constants [L ,x0, k, b]
    '''
    time_series = exp_dict[tf][gene]
    xdata = []
    ydata = []
    for time in time_series:
        xdata.append(time[0])
        ydata.append(time[1])
    ydata = scale_data(ydata)
    p0 = [max(ydata), np.median(xdata), 1, min(ydata)]    
    
    try:
        # Curve fitting with bounds and method specified
        bounds = ([-np.inf, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf])
        params, _ = curve_fit(sigmoid, xdata, ydata, p0=p0, bounds=bounds, method='dogbox', maxfev=100000)
        return params
    # except Exception:
    #     pass
    # try:
    #     bounds = ([-np.inf, 0, 0, -np.inf, -np.inf, 0, 0, -np.inf], 
    #               [np.inf, np.inf, np.inf, np.inf, -np.inf, 0, 0, -np.inf])
    #     params, _ = curve_fit(double_sigmoid, xdata, ydata, p0=p0, bounds=bounds, method='dogbox', maxfev=100000)
    #     return params
    except Exception as e:
        print(f"{e} - TF: {tf}, target: {gene}")
        return None

# building heat maps
def heat_map_colors():
    '''
    Set the color scheme of heat map to Haase lab colors.
    
    Returns:
    haase (matplotlib.colors.LinearSegmentedColormap): color scheme
    '''
    norm = matplotlib.colors.Normalize(-1.5,1.5)
    colors = [[norm(-1.5), "cyan"],
          [norm(0), "black"],
         [norm(1.5), "yellow"]]
    haase = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    return haase

def create_heat_maps(df):
    '''
    Build heat maps based on gene expression levels, and creates .png file containing output. 
    Arrange subplots based on TF. 

    Parameters:
    df (pandas DataFrame): gene expression data of desired TFs and targets
    '''
    # run heat_map_colors() to get haase color scheme
    haase = heat_map_colors()

    # num_cols and num_rows to determine dimensions of plot and subplots
    tfs = df['TF'].unique()
    num_tfs = len(tfs)
    num_cols = 4
    num_rows = math.ceil(num_tfs / num_cols)

    # create plot
    fig = plt.figure(figsize = (15,10 + num_rows * 3))
    fig.subplots_adjust(hspace=0.4, wspace=0.4, top = 0.90)
    fig.suptitle("IDEA Dataset Expression Levels", fontsize = 15)
    

    for i, tf in enumerate(tfs):    
        # Filter data for the specific TF
        tf_data = df[df['TF'] == tf]

        # Create a pivot table for heatmap 
        heatmap_data = tf_data.pivot(index='GeneName', columns='time', values='log2_cleaned_ratio')

        # Make subplot
        ax = plt.subplot(num_rows, num_cols, i + 1)
        sns.heatmap(heatmap_data, cmap=haase, cbar=True, vmin=-2, vmax=2, ax=ax)
        ax.set_title(tf)
        ax.set_xlabel('time (min)')
        ax.set_ylabel('')  

        # set tick positions and labels
        # times = np.array([0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 90.0])
        gene_names = heatmap_data.index.to_numpy()
        # ax.set_xticks(np.arange(len(times)), labels=times)
        ax.set_yticks(np.arange(len(gene_names)) + 0.25, labels=gene_names, fontsize=6)
        plt.xticks(rotation=45)  

    plt.savefig("heat_maps.png")
    plt.show()

# network construction
def build_ref_network():
    df = pd.read_csv('ref_edges.csv')
    for _, row in df.iterrows():
        tf = row['reg']
        if tf not in gene_nodes:
            gene_nodes[tf] = structs.GeneNode(tf)
        tf_node = gene_nodes[tf]

        target = row['target']
        if target not in gene_nodes:
            gene_nodes[target] = structs.GeneNode(target)
        target_node = gene_nodes[target]

        act = row['type'] == 'act'

        edge = tf_node.get_edge(target_node, act)
        if edge is None:
            tf_node.add_edge(target_node, act)


def build_tree(tf, t_thresh, amp_thresh, edges_df):
    '''
    Create a tree structure with a height of 1 for a given TF.

    Parameters:
    tf (string): TF
    t_thresh (float): time threshold for direct connection
    amp_thresh (int or float): minimum amplitude threshold for activation/inhibition
    edges_df (pandas DataFrame): running DataFrame of current Edges
    '''
    # initializing the root node if it doesn't exist already
    if tf not in gene_nodes:
        gene_nodes[tf] = structs.GeneNode(tf)
    
    node = gene_nodes[tf]

    # establish edges  
    for gene in exp_dict[tf].keys(): # target
        time, sign = get_sig_info(tf, gene, amp_thresh)
        if time is not None and time > 0 and time < t_thresh:
            if gene not in gene_nodes:
                gene_nodes[gene] = structs.GeneNode(gene)
            target = gene_nodes[gene]
            node.add_edge(target, sign) 
            edges_df.loc[len(edges_df.index)] = [tf, gene, 'act' if sign else 'rep']
        
def visualize_gene_network(timestamp):
    '''
    Create visualization of GRN using graphviz package. 
    Output stored in a .pdf file.
    '''
    # produce ref graph
    build_ref_network()

    # create digraph
    dot = graphviz.Digraph(comment='Gene Regulatory Network')
    
    # Add nodes
    for gene in gene_nodes:
        dot.node(gene_nodes[gene].gene, shape='box')
    
    # Add edges
    both, refs, grns = compare_edges('ref_edges.csv', 'grn_edges.csv')
    for gene in gene_nodes:
        for edge in gene_nodes[gene].edges:
            # determine arrowhead
            if edge.act:
                arrowhead = 'normal'
            else:
                arrowhead = 'tee'
            # determine edge color
            if edge in both:
                color = '#008000'
            elif edge in refs:
                color = '#0059b3'
            elif edge in grns:
                color = '#b30000'
            else:
                print("edge error")
            dot.edge(gene, edge.target.gene, color=color, arrowhead=arrowhead)

    # Render the graph
    dot.render(f'grns/gene_network_{timestamp}', view=True)

def build_network(df):
    '''
    The primary function for building the overall network.

    Parameters:
    df (pandas DataFrame): contains gene expression data of desired TFs and targets
    '''
    # get timestamp
    ct = datetime.datetime.now()
    timestamp = ct.strftime("%Y-%m-%d %H:%M:%S")

    # build network
    extract_expression_data(df)
    tfs = df['TF'].unique()
    edges_df = pd.DataFrame(columns=['reg', 'target', 'type'])
    thresh_df = pd.DataFrame(columns=['reg', 't_thresh', 'amp_thresh'])
    for tf in tfs:
        t_thresh = 15
        amp_thresh = 0.5
        build_tree(tf, t_thresh, amp_thresh, edges_df) # dummy thresholds
        thresh_df.loc[len(thresh_df.index)] = [tf, t_thresh, amp_thresh]
    edges_df.to_csv("grn_edges.csv", index=False)
    thresh_df.to_csv(f"grn_params/params_{timestamp}.csv", index=False)

    visualize_gene_network(timestamp)

# compare with reference graph
def compare_edges(ref, grn):
    '''
    Compare edges in ref_edges.csv and grn_edges.csv to find agreements and discrepancies.

    Parameters:
    ref (string): path to ref_edges.csv
    grn (string): path to grn_edges.csv

    Returns: 
    both (list of Edges): Edges found in both files
    refs (list of Edges): Edges found only in ref_edges.csv and not in grn_edges.csv
    grns (list of Edges): Edges found only in grn_edges.csv and not in ref_edges.csv
    '''
    # read csv files into DataFrames
    ref_df = pd.read_csv(ref)
    grn_df = pd.read_csv(grn)

    # get common edges
    common_rows = pd.merge(ref_df, grn_df, how='inner')
    both = df_to_edges(common_rows)

    # get edges that are diff
    ref_rows = pd.merge(ref_df, grn_df, how='left').drop_duplicates(keep=False)
    refs = df_to_edges(ref_rows)

    grn_rows = pd.merge(ref_df, grn_df, how='right').drop_duplicates(keep=False)
    grns = df_to_edges(grn_rows)

    # return lists of Edges
    return both, refs, grns

def df_to_edges(df):
    '''
    Given a pandas DataFrame containing edge information, return the Edges listed.

    Parameters: 
    df (pandas DataFrame): input DataFrame with columns 'reg', 'target', and 'type'

    Returns:
    edges (list of Edges): the list of Edges corresponding to the input df
    '''
    edges = []
    for _, row in df.iterrows():
        reg = gene_nodes[row['reg']] if row['reg'] in gene_nodes else None
        target = gene_nodes[row['target']] if row['target'] in gene_nodes else None
        if reg == None or target == None:
            break
        act = row['type'] == 'act'
        edge = reg.get_edge(target, act)
        edges.append(edge)
    return edges

# main function
def main():
    dir = "/Users/jingliu/Documents/haase/IDEA_data"
    file = "idea_tall_expression_data.tsv"
    path = os.path.join(dir, file)

    gene_dir = '/Users/jingliu/Documents/haase/haase_lab'
    gene_file = 'genes_sorted.csv'
    gene_path = os.path.join(gene_dir, gene_file)

    df = filter_df(path, gene_path)

    # # create heat maps
    # create_heat_maps(df)

    # reformat data and build network
    build_network(df)

if __name__ == '__main__':
    main()   