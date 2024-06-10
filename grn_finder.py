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

importlib.reload(structs)

# global list of GeneNode names that already exist
existing_nodes = []

# global dictionary of expression data
exp_dict = {}
# 
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
        exp_dict[tf][gene].append((row['time'], row['log2_cleaned_ratio']))

def get_t_act(row):
    '''
    row (pandas Series): row from the original DataFrame
    return: time of activation/inhibition, and whether its activation or not (boolean)
    '''
    tf = row['TF']
    gene = row['GeneName']
    time_series = exp_dict[tf][gene]
    opt = sigmoid_curve_fit(time_series) # [L_opt, x0_opt, k_opt, b_opt]
    return (opt[0]/2.0 + opt[3]), (opt[0] > 0)

# borrowed from online: https://stackoverflow.com/questions/55725139/fit-sigmoid-function-s-shape-curve-to-data-using-python
def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return y

def sigmoid_curve_fit(time_series):
    '''
    time_series: list of data from log2_cleaned_ratio
    '''
    xdata = [0, 5, 10, 15, 20, 30, 45, 90]
    ydata = time_series
    p0 = [max(ydata), np.median(xdata), 1, min(ydata)] # this is an mandatory initial guess
    opt, cov = scipy.optimize.curve_fit(sigmoid, xdata, ydata,p0, method='dogbox')
    return opt

def build_tree(tf, df, threshold):
    '''
    Create a tree structure with a height of 1 for a given TF.
    tf (string): TF
    df (DataFrame): Pandas DataFrame with expression data
    '''
    node_name = "node_" + tf

    # initializing the root node if it doesn't exist already
    if tf not in existing_nodes:
        df_tf = df[df['TF'] == tf]
        globals()[node_name] = structs.GeneNode(tf)
        existing_nodes.append(globals()[node_name])
    
    node = globals()[node_name]

    # establish edges
    for _, row in df_tf.iterrows():
        # falls = get_t_fall(row)
        # rises = get_t_rise(row)
        target = row['GeneNode']
        # time = max(falls[0], rises[0])
        time, act = get_t_act(row)
        if time < threshold:
            edge = structs.Edge(target, act)
            node.add_edge(edge)    
        
def visualize_gene_network(gene_nodes):
    dot = graphviz.Digraph(comment='Gene Regulatory Network')
    
    # Add nodes
    for gene in gene_nodes:
        dot.node(gene.gene)
    
    # Add edges
    for gene in gene_nodes:
        for edge in gene.edges:
            if edge.act:
                color = 'black'
                arrowhead = 'normal'
            else:
                color = 'black'
                arrowhead = 'tee'
            dot.edge(gene.gene, edge.target.gene, color=color, arrowhead = arrowhead)
    
    # Render the graph
    dot.render('gene_network', view=True)

def build_network(tfs, df):
    '''
    The primary function for building the overall network.
    tfs (list): list of strings representing TFs
    '''
    for tf in tfs:
        build_tree(tf, df, 20) # dummy threshold

    visualize_gene_network(existing_nodes)

def main():
    dir = "/Users/jingliu/Documents/haase/IDEA_data"
    file = "idea_tall_expression_data.tsv"
    path = os.path.join(dir, file)
    df = pd.read_csv(path, sep='\t')
    tfs = ['ACA1', 'ACE2']
    filtered_df = df[df['TF'].isin(tfs)]
    extract_expression_data(filtered_df)
    build_network(tfs, df)

if __name__ == '__main__':
    main()
    
    