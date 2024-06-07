import pandas as pd
import os
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import glob
import numpy as np

def extract_expression_data(path, sep=',', strain=None):
    '''
    This function extracts the necessary information from 
    path(string): includes location and name of file
    sep (string): separator in file (e.g. "\t" for .tsv files, "," for .csv files)
    return: 
    '''
    df = pd.read_csv(path, sep=sep)
    if strain is not None: 
        df = df[df['strain'] == strain]
    df = df[['TF', 'strain', 'GeneName', 'time', 'log2_shrunken_timecourses']]

def extract_TF_data(df, tf):
    '''
    df (Pandas DataFrame): DataFrame containing the expression data of all TFs
    tf (string)
    '''

def get_t_fall():
    '''
    return: list of floats
    '''
    
def get_t_rise():
    '''
    return: list of floats
    '''

def build_tree():
    '''
    '''