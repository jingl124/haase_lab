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
import datetime
import ast
import grn_search as search
import grn_finder as finder

importlib.reload(structs)
importlib.reload(finder)

groups = {}
# key: group ID
# value: GroupNode

# Group node creation and interactions
def create_group_nodes(orthologs, complexes):
    """
    Creates group nodes (orthologs and complexes) and adds edges between them based on the specified criteria.

    Args:
        orthologs (list): List of ortholog group nodes.
        complexes (list): List of complex group nodes.

    Returns:
        None
    """

    for group in orthologs:
        for node in group.nodes:
            for edge in node.edges:
                if group.get_edge(edge.target, edge.act) is None:
                    group.add_edge(edge.target, edge.act)

    for group in complexes:
        for node in group.nodes:
            for edge in node.edges:
                if group.get_edge(edge.target, edge.act) is None:
                    group.add_edge(edge.target, edge.act)

    # Update gene_nodes with group nodes
    for group in orthologs + complexes:
        groups[group.name] = group

def in_ortholog(group):
    '''
    For a given group of orthologs, compile all the in edges based on the nodes in the group.
    '''
    list = []
    for node in group.nodes:
        starts = search.get_in_edges(node.name)
        for s in starts: # s is tuple (target.name, act)
            start = group[s[0]]
            if start.get_edge(group, s[1]) == None:
                start.add_edge(group, s[1])
                list.append((start, s[1])) # starting node and activation
    return list

def read_group_nodes():
    '''
    Returns:
    orthologs (list of GroupNodes):
    complexes (list of GroupNodes):

    '''
    df = pd.read_csv("ref_groups.csv")

    orthologs = []
    complexes = []
    for _, row in df.iterrows():
        group = row['group']
        type = row['type']
        nodes_str = group[1:len(group)-1].split(",")
        nodes = []
        for str in nodes_str:
            nodes.append(gene_nodes[str])
        node = structs.GroupNode(nodes, type)
        if type == 'ortholog':
            orthologs.append(node)
        elif type == 'complex':
            complexes.append(node)
        gene_nodes[node.name] = node
    
    return orthologs, complexes

def main():
    orthologs, complexes = 

if __name__ == '__main__':
    main() 
