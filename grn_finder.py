import pandas as pd
import os
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

# global list of Edges that already exist
# 

def extract_expression_data(path, sep=','):
    '''
    This function extracts the necessary information from 
    path(string): includes location and name of file
    sep (string): separator in file (e.g. "\t" for .tsv files, "," for .csv files)
    return: 
    '''
    df = pd.read_csv(path, sep=sep)
    # selected_columns = ['TF', 'strain', 'GeneName', 'time', 'log2_shrunken_timecourses']
    # df = df[selected_columns]
    print("Expression data extracted.")
    return df

def get_t_fall(row):
    '''
    row (pandas Series): row from the original DataFrame
    return: list of floats
    '''

    
def get_t_rise(row):
    '''
    return: list of floats
    '''

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
        falls = get_t_fall(row)
        rises = get_t_rise(row)
        target = row['GeneNode']
        time = max(falls[0], rises[0])
        act = (time == rises[0])
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

def build_network(tfs):
    '''
    The primary function for building the overall network.
    tfs (list): list of strings representing TFs
    '''
    for tf in tfs:
        build_tree(tf)

    visualize_gene_network(existing_nodes)
    
    