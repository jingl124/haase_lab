# functions to search within the built grn and paths
import pandas as pd
import os
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import numpy as np
import grn_finder as grn

def path_search(start=None, end=None, sign=None):
    '''
    Return strings of paths that satisfy the given arguments.

    Parameters:
    start (string): 
    end (string):
    sign (string): 'act' if activating, 'rep' if repressing. 
        Returns both activating and repressing edges if sign == None.

    
    '''
    # invalid arguments
    if start == None and end == None:
        raise ValueError("Needs a start or end node.")
    
    # read in .csv file
    path_df = pd.read_csv("paths.csv")
    filtered_df = path_df.copy()

    # filtering based on arguments
    if start is not None:
        start = start.upper()
        filtered_df = filtered_df[filtered_df['start'] == start]
    if end is not None:
        end = end.upper()
        filtered_df = filtered_df[filtered_df['end'] == end]
    if sign is not None:
        filtered_df = filtered_df[filtered_df['sign'] == sign]

    # no rows left
    if len(filtered_df.index) == 0:
        print("No paths that satisfy arguments.")
        return

    # print out paths
    paths = ""
    for _, row in filtered_df.iterrows():
        path = row['path']
        paths += f"{path}\n"
    print(paths)

def get_out_edges(node_name):
    '''
    Given the name of a gene, find all edges pointing from this node.

    Parameters: 
    node_name (string): name of the node

    Returns:
    nodes (list of strings): target names of all the edges
    '''
    node = grn.gene_nodes[node_name]
    nodes = []
    for edge in node.edges:
        nodes.append(edge.target.gene)
    print(nodes)
    return nodes

def get_in_edges(node_name):
    '''
    Given the name of a gene, find all edges that point to this node.

    Parameters: 
    node_name (string): name of the node

    Returns:
    nodes (list of strings): target names of all the edges
    '''
    target = grn.gene_nodes[node_name]
    nodes = []
    for gene in grn.gene_nodes.keys():
        node = grn.gene_nodes[gene]
        edge = node.get_edge(target)
        if edge == None:
            continue
        else:
            nodes.append(node.gene)
    print(nodes)
    return nodes

def search_paths(file, start, end):
    '''
    Given a file containing paths generated from grn_finder.py, search for paths with the specified start and end nodes.

    Parameters:
    file (string): file path
    start (string): start node
    end (string): end node

    Returns:
    df (pandas DataFrame): DataFrame including all paths satisfying the specified start and end nodes
    '''
    df = pd.read_csv(file)



def main():
    '''
    '''
    grn.main()

if __name__ == '__main__':
    main()

