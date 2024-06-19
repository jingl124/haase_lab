import pandas as pd
import os
import scipy.optimize
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

def construct_nodes():
    '''
    Constructing the reference cell cycle network into GeneNodes via hard code.
    GeneNodes constructed will be accessible in the global dictionary gene_nodes. 
    '''

def build_graph():
    '''
    Build a graphviz graph based on the constructed GeneNodes.
    '''