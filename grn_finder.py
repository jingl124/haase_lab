import pandas as pd
import os
import scipy.optimize
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import glob
import numpy as np
import structs
import importlib
import graphviz

importlib.reload(structs)

# global list of GeneNode names (strings) that already exist
existing_nodes = []

# global dictionary of expression data
exp_dict = {}
# 
def extract_expression_data(df):
    '''
    This function extracts the necessary information from 
    df (pandas DataFrame): contains all data
    '''
    for _, row in df.iterrows():
        tf = row['TF']
        gene = row['GeneName']
        if tf not in exp_dict:
            exp_dict[tf] = {}  # Initialize nested dictionary for tf if not present
        if gene not in exp_dict[tf]:
            exp_dict[tf][gene] = []  # Initialize list for gene if not present
        exp_dict[tf][gene].append(row['log2_cleaned_ratio'])

def get_t_act(tf, gene, amp_thresh):
    '''
    row (pandas Series): row from the original DataFrame
    return: time of activation/inhibition, and whether its activation or not (boolean)
    '''
    time_series = exp_dict[tf][gene]
    opt = sigmoid_curve_fit(time_series) # [L_opt, x0_opt, k_opt, b_opt]
    if opt is None or len(opt) == 0 or abs(opt[0]) < amp_thresh:
        return None
    return [opt[1], (opt[0] > 0)]

# borrowed from online: https://stackoverflow.com/questions/55725139/fit-sigmoid-function-s-shape-curve-to-data-using-python
def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return y

def scale_data(data):
    """
    Scale the data to a range suitable for fitting.
    """
    data_min = min(data)
    data_max = max(data)
    if data_max == data_min:
        return [0 for _ in data]
    return [(d - data_min) / (data_max - data_min) for d in data]

def sigmoid_curve_fit(time_series):
    '''
    time_series: list of data from log2_cleaned_ratio
    '''
    xdata = [0, 5, 10, 15, 20, 30, 45, 90]
    ydata = scale_data(time_series)
    # p0 = [0, 0, 0, 0] #[max(ydata), np.median(xdata), 1, min(ydata)] # this is an mandatory initial guess
    # opt, cov = scipy.optimize.curve_fit(sigmoid, xdata, ydata, p0, method='dogbox', maxfev=100000)
    #  # Initial guess for the parameters
    # p0 = [max(ydata), np.median(xdata), 1, min(ydata)]
    
    # # Bounds for the parameters to ensure all x and y values are > 0
    # bounds = (0, [np.inf, np.inf, np.inf, np.inf])
    
    # # Curve fitting with bounds and method specified
    # opt, cov = scipy.optimize.curve_fit(sigmoid, xdata, ydata, p0=p0, bounds=bounds, method='dogbox', maxfev=100000)
    # Initial guess for the parameters
    p0 = [max(ydata), np.median(xdata), 1, min(ydata)]
    
    # Bounds for the parameters to ensure positive x0, and k, but allow L and b to be any value
    bounds = ([-np.inf, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf])
    
    try:
        # Curve fitting with bounds and method specified
        opt, _ = scipy.optimize.curve_fit(sigmoid, xdata, ydata, p0=p0, bounds=bounds, method='dogbox', maxfev=100000)
        return opt
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
    # return opt

def build_tree(tf, df, t_thresh, amp_thresh):
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
        existing_nodes.append(node_name)
    
    node = globals()[node_name]

    # establish edges
    for gene in exp_dict[tf].keys():#_, row in df_tf.iterrows():
        # falls = get_t_fall(row)
        # rises = get_t_rise(row)
        # time = max(falls[0], rises[0])
        time = get_t_act(tf, gene, amp_thresh)
        if time is not None and time[0] > 0 and time[0] < t_thresh:
            target_name = "node_" + gene
            if gene not in existing_nodes:
                globals()[target_name] = structs.GeneNode(gene)
                existing_nodes.append(target_name)
            target = globals()[target_name]
            node.add_edge(target, time[1])    
        
def visualize_gene_network(gene_nodes):
    dot = graphviz.Digraph(comment='Gene Regulatory Network')
    
    # Add nodes
    for gene in gene_nodes:
        dot.node(globals()[gene].gene)
    
    # Add edges
    for gene in gene_nodes:
        for edge in globals()[gene].edges:
            if edge.act:
                color = 'black'
                arrowhead = 'normal'
            else:
                color = 'black'
                arrowhead = 'tee'
            dot.edge(globals()[gene].gene, edge.target.gene, color=color, arrowhead = arrowhead)
    
    # Render the graph
    dot.render('gene_network', view=True)
    # G = nx.DiGraph()
    
    # # Add nodes and edges
    # for gene_name in gene_nodes:
    #     gene_node = globals()[gene_name]
    #     G.add_node(gene_node.gene)
    #     for edge in gene_node.edges:
    #         if edge.act:
    #             color = 'green'
    #             arrowhead = 'normal'
    #         else:
    #             color = 'red'
    #             arrowhead = 'tee'
    #         G.add_edge(gene_node.gene, edge.target.gene, color=color, arrowhead=arrowhead)
    
    # # Get edge colors and styles
    # edge_colors = [G[u][v]['color'] for u, v in G.edges()]
    # edge_styles = ['solid' if G[u][v]['arrowhead'] == 'normal' else 'dashed' for u, v in G.edges()]

    # # Draw the network
    # pos = nx.spring_layout(G)  # Position nodes using Fruchterman-Reingold force-directed algorithm
    
    # nx.draw(G, pos, with_labels=True, node_size=3000, node_color='lightblue', font_size=10, font_weight='bold', edge_color=edge_colors, style=edge_styles, arrowsize=20)
    # plt.title('Gene Regulatory Network')
    # plt.show()

def build_network(tfs, df):
    '''
    The primary function for building the overall network.
    tfs (list): list of strings representing TFs
    '''
    for tf in tfs:
        build_tree(tf, df, 20, 0.2) # dummy thresholds

    visualize_gene_network(existing_nodes)

def main():
    dir = "/Users/jingliu/Documents/haase/IDEA_data"
    file = "idea_tall_expression_data.tsv"
    path = os.path.join(dir, file)
    df = pd.read_csv(path, sep='\t')
    tfs = ['ACA1']
    filtered_df = df[df['TF'].isin(tfs)]
    filtered_df = filtered_df.iloc[:180]
    extract_expression_data(filtered_df)
    build_network(tfs, df)

if __name__ == '__main__':
    main()
    
    