import pandas as pd
import os
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import glob
import numpy as np

def dummy():
    print("hi")

def extract_expression_data(path, sep=','):
    '''
    This function extracts the necessary information from 
    path(string): includes location and name of file
    sep (string): separator in file (e.g. "\t" for .tsv files, "," for .csv files)
    return: 
    '''
    df = pd.read_csv(path, sep=sep)
    selected_columns = ['TF', 'strain', 'GeneName', 'time', 'log2_shrunken_timecourses']
    df = df[selected_columns]
    print("Expression data extracted.")
    return df

def get_t_fall():
    '''
    return: list of floats
    '''
    
def get_t_rise():
    '''
    return: list of floats
    '''

def build_tree(tf, df, threshold):
    '''
    Create a tree structure with a height of 1 for a given TF.
    tf (string): TF
    df (DataFrame): Pandas DataFrame with expression data
    '''
    # df_tf = df[]

def build_network(tfs):
    '''
    The main function for building the overall network.
    tfs (list): list of strings representing TFs
    '''
    for tf in tfs:
        build_tree(tf)
    