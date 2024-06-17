import pandas as pd
import os
import scipy.optimize
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import glob
import numpy as np
import structs
import importlib
import graphviz
import matplotlib
import math

importlib.reload(structs)

# global dictionary of GeneNode names (strings) and the actual GeneNodes that already exist
gene_nodes = {}

# global dictionary of expression data
exp_dict = {}

# extracting and wrangling data
def extract_expression_data(df):
    '''
    This function extracts the necessary information from 
    df (pandas DataFrame): contains all data
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

# using sigmoidal fit to retrieve gene expression info
def get_t_act(tf, gene, amp_thresh):
    '''
    row (pandas Series): row from the original DataFrame
    return: time of activation/inhibition, and whether its activation or not (boolean)
    '''
    opt = sigmoid_curve_fit(tf, gene) # [L_opt, x0_opt, k_opt, b_opt]
    if opt is None or len(opt) == 0 or abs(opt[0]) < amp_thresh:
        return None, None
    return opt[1], (opt[0] > 0)

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return y

def scale_data(data):
    """
    Scale the data to a range suitable for fitting.
    """
    data_min = min(data)
    data_max = max(data)
    if data_max == data_min:
        return [0 for _ in data]
    return [(d - data_min) / (data_max - data_min) for d in data]

def sigmoid_curve_fit(tf, gene):
    '''
    time_series: list of data from log2_cleaned_ratio
    '''
    time_series = exp_dict[tf][gene]
    xdata = []
    ydata = []
    for time in time_series:
        xdata.append(time[0])
        ydata.append(time[1])
    ydata = scale_data(ydata)
    p0 = [max(ydata), np.median(xdata), 1, min(ydata)]
    
    # Bounds for the parameters to ensure positive x0, and k, but allow L and b to be any value
    bounds = ([-np.inf, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf])
    
    try:
        # Curve fitting with bounds and method specified
        opt, _ = scipy.optimize.curve_fit(sigmoid, xdata, ydata, p0=p0, bounds=bounds, method='dogbox', maxfev=100000)
        return opt
    except Exception as e:
        print(f"An error occurred - TF: {tf}, target: {gene}")
        return None

# building heat maps
def heat_map_colors():

    norm = matplotlib.colors.Normalize(-1.5,1.5)
    colors = [[norm(-1.5), "cyan"],
          [norm(0), "black"],
         [norm(1.5), "yellow"]]
    haase = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    return haase

def create_heat_maps(df):
    # run heat_map_colors() to get haase color scheme
    haase = heat_map_colors()

    # create plot
    fig = plt.figure(figsize = (15,10))
    fig.subplots_adjust(hspace=0.4, wspace=0.4, top = 0.90)
    fig.suptitle("IDEA Dataset Expression Levels", fontsize = 15)

    tfs = df['TF'].unique()
    num_tfs = len(tfs)

    # Calculate number of rows needed for subplots (2 columns per row)
    num_cols = 4
    num_rows = math.ceil(num_tfs / num_cols)

    for i, tf in enumerate(tfs):    
        # Filter data for the specific TF
        tf_data = df[df['TF'] == tf]

        # Create a pivot table for heatmap 
        heatmap_data = tf_data.pivot(index='GeneName', columns='time', values='log2_cleaned_ratio')

        # Make subplot
        plt.subplot(num_rows, num_cols, i + 1)
        sns.heatmap(heatmap_data, cmap=haase, cbar=True)
        plt.title(tf)

    plt.savefig("heat_maps.png")
    plt.show()




# network construction
def build_tree(tf, t_thresh, amp_thresh):
    '''
    Create a tree structure with a height of 1 for a given TF.
    tf (string): TF
    t_thresh: time threshold for direct connection
    amp_thresh: minimum amplitude threshold for activation/inhibition
    '''
    # initializing the root node if it doesn't exist already
    if tf not in gene_nodes:
        gene_nodes[tf] = structs.GeneNode(tf)
    
    node = gene_nodes[tf]

    # establish edges  
    for gene in exp_dict[tf].keys():
        time, sign = get_t_act(tf, gene, amp_thresh)
        if time is not None and time > 0 and time < t_thresh:
            if gene not in gene_nodes:
                gene_nodes[gene] = structs.GeneNode(gene)
            target = gene_nodes[gene]
            node.add_edge(target, sign) 
        
def visualize_gene_network(gene_nodes):
    dot = graphviz.Digraph(comment='Gene Regulatory Network')
    
    # Add nodes
    for gene in gene_nodes:
        dot.node(gene_nodes[gene].gene, shape='box')
    
    # Add edges
    edges_df = pd.DataFrame(columns=['TF', 'GeneName', 'arrowhead'])
    for gene in gene_nodes:
        for edge in gene_nodes[gene].edges:
            if edge.act:
                arrowhead = 'normal'
            else:
                arrowhead = 'tee'
            dot.edge(gene, edge.target.gene, color='black', arrowhead=arrowhead)
    edges_df.to_csv("edges.csv", index=False)

    # Render the graph
    dot.render('gene_network', view=True)

def build_network(df):
    '''
    The primary function for building the overall network.
    tfs (list): list of strings representing TFs
    '''
    tfs = df['TF'].unique()
    for tf in tfs:
        build_tree(tf, 20, 0.2) # dummy thresholds

    visualize_gene_network(gene_nodes)

def main():
    dir = "/Users/jingliu/Documents/haase/IDEA_data"
    file = "idea_tall_expression_data.tsv"
    path = os.path.join(dir, file)

    gene_dir = '/Users/jingliu/Documents/haase/haase_lab'
    gene_file = 'genes_sorted.csv'
    gene_path = os.path.join(gene_dir, gene_file)

    df = filter_df(path, gene_path)

    # # reformat data and build network
    # extract_expression_data(df)
    # build_network(df)

    # create heat maps
    create_heat_maps(df)

if __name__ == '__main__':
    main()
    
    