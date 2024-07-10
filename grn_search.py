# functions to search within the built grn and paths
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

path_search(start='ACE2', end='ACE2')