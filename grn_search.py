# functions to search within the built grn and paths
import pandas as pd
import os
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import numpy as np
import grn_finder as grn
import datetime
import sys
import argparse

# global timestamp
timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

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

def search_paths(file, start=None, end=None):
    '''
    Given a pandas DataFrame containing all paths found in a GRN, search for paths with the specified start and end nodes.

    Parameters:
    df (pandas DataFrame): contains all the paths
    start (string): start node
    end (string): end node

    Returns:
    paths_df (pandas DataFrame): DataFrame including all paths satisfying the specified start and end nodes
    '''
    df = pd.read_csv(file)

    if 'start' not in df.columns:
        raise ValueError("File does not contain 'start' column; try again.")
    elif 'end' not in df.columns:
        raise ValueError("File does not contain 'end' column; try again.")
    
    paths_df = df[(df['start'] == start) & (df['end'] == end)]
    
    if paths_df.empty:
        print("No paths found for the given start and end points.")
    else:
        print("Paths found:")
        print(paths_df)

    return paths_df

def find_ffls(df):
    '''
    Given a pandas DataFrame containing all paths found in a GRN, find all FFLs and return them in a string. 

    Parameters: 
    df (pandas DataFrame): DataFrame containing all paths

    Returns:
    str (string): string to be put into an output text file
    '''
    str = ""

    starts = df['start'].unique()
    ends = df['end'].unique()

    for start in starts:
        for end in ends:
            paths = df[(df['start'] == start) & (df['end'] == end)]
            if paths.empty or len(paths) <= 1:
                continue
            str += f"{start}, {end}:\n{paths.to_string(index=False)}\n"

    return str

def sample_function():
    print("hello world")

def main():
    # grn.main()

    parser = argparse.ArgumentParser(description='GRN search utility')
    parser.add_argument('function', type=str, help='Function name to execute')
    parser.add_argument('--args', nargs='*', help='Arguments for the function', default=[])

    args = parser.parse_args()
    function_name = args.function
    function_args = args.args

    # Check if the function exists in the current module
    if function_name in globals():
        func = globals()[function_name]
        if callable(func):
            func(*function_args)
        else:
            print(f"{function_name} is not a callable function.")
    else:
        print(f"Function {function_name} not found.")

if __name__ == '__main__':
    main()

