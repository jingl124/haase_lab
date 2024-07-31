# functions to search within the built grn and paths
import pandas as pd
import os
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import numpy as np
import grn_finder as grn
import datetime
import argparse

# global timestamp
timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

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
            nodes.append((node.gene, edge.act))
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
    
    paths_df = df.copy()
    
    if start != 'None' and start != None:
        paths_df = paths_df[paths_df['start'] == start]
    if end != 'None' and end != None:
        paths_df = paths_df[paths_df['end'] == end]
    
    if paths_df.empty:
        print("No paths found.")
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
    grn.main()

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
            # Parse key-value arguments or positional arguments
            parsed_args = {}
            positional_args = []

            for arg in function_args:
                if '=' in arg:
                    key, value = arg.split('=', 1)
                    parsed_args[key] = value
                else:
                    positional_args.append(arg)

            # Check for None values indicated by empty strings
            parsed_args = {k: None if v == '' else v for k, v in parsed_args.items()}

            # Call the function with parsed arguments
            if parsed_args:
                # If there are keyword arguments, use them
                func(**parsed_args)
            else:
                # Use positional arguments if no keyword arguments are found
                func(*positional_args)
        else:
            print(f"{function_name} is not a callable function.")
    else:
        print(f"Function {function_name} not found.")

if __name__ == '__main__':
    main()

