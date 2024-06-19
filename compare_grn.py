import pandas as pd
import os
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import numpy as np
import structs
import importlib
import graphviz
import matplotlib
import math

importlib.reload(structs)

# global dictionary of GeneNode names (strings) and the actual GeneNodes that already exist
gene_nodes = {}

def construct_nodes(path):
    '''
    Create GeneNodes from a list of genes, read in from input .txt file. 
    GeneNodes constructed will be accessible in the global dictionary gene_nodes. 

    Parameters:
    path (string): location of input .txt file containing genes
    '''
    file = open(path, 'r')
    nodes = file.readlines()
    for node in nodes:
        node = node.strip()
        gene_nodes[node] = structs.GeneNode(node)
    
def construct_edges(path, sep=','):
    '''
    Create Edges from a list of edges, read in from input .csv file.

    Parameters:
    path (string): location of input .csv file containing edges
    sep (string): separator for reading in input file, default is for .csv
    '''
    edges_df = pd.read_csv(path, sep=sep)
    for _, row in edges_df.iterrows():
        reg = gene_nodes[row['reg']]
        target = gene_nodes[row['target']]
        act = row['type'] == 'act'
        reg.add_edge(target, act)


def build_ref_graph():
    '''
    Build a graphviz graph based on constructed GeneNodes and Edges.
    '''
    # build GeneNode data structure
    dir = "/Users/jingliu/Documents/haase/haase_lab"
    node_path = os.path.join("Cell_cycle_GRN_gene.txt")
    edge_path = os.path.join(dir, "ref_edges.csv")
    construct_nodes(node_path)
    construct_edges(edge_path)

    # creating graphviz graph
    dot = graphviz.Digraph(comment='Cell Cycle Reference Network')
    
    # Add nodes
    for gene in gene_nodes:
        dot.node(gene_nodes[gene].gene, shape='box')
    
    # Add edges
    for gene in gene_nodes:
        for edge in gene_nodes[gene].edges:
            if edge.act:
                arrowhead = 'normal'
            else:
                arrowhead = 'tee'
            dot.edge(gene, edge.target.gene, color='black', arrowhead=arrowhead)

    # Render the graph
    dot.render('ref_network', view=True)

def main():
    build_ref_graph()

if __name__ == '__main__':
    main()