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

importlib.reload(structs)

# global dictionary of GeneNode names (strings) and the actual GeneNodes that already exist
gene_nodes = {}

# global dictionary of expression data
exp_dict = {}

# global pandas DataFrame for time and amplitude thresholds
thresh_df = pd.DataFrame(columns=['gene', 't_thresh', 'amp_thresh'])

# global timestamp
timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# Constants
DATA_DIR = "/Users/jingliu/Documents/haase/IDEA_data"
DATA_FILE = "idea_tall_expression_data.tsv"
GENE_DIR = '/Users/jingliu/Documents/haase/haase_lab'
GENE_FILE = 'genes_info.csv'


# extracting and wrangling data
def extract_expression_data(df):
    '''
    This function extracts the necessary information from the DataFrame and stores it into global dictionary exp_dict.
    exp_dict is formatted like this:
        {
            TF1: {
                target_gene1: [(time, log2_cleaned_ratio), (time, log2_cleaned_ratio), ...],
                target_gene2: [(time, log2_cleaned_ratio), (time, log2_cleaned_ratio), ...],
                ...
            }
            TF2: {
                target_gene1: [(time, log2_cleaned_ratio), (time, log2_cleaned_ratio), ...],
                target_gene2: [(time, log2_cleaned_ratio), (time, log2_cleaned_ratio), ...],
                ...
            }
        }
    Parameters:
    df (pandas DataFrame): contains all data derived from gene expression data file
    '''
    for _, row in df.iterrows():
        tf = row['TF']
        gene = row['GeneName']
        if tf not in exp_dict:
            exp_dict[tf] = {}  # Initialize nested dictionary for tf if not present
        if gene not in exp_dict[tf]:
            exp_dict[tf][gene] = []  # Initialize list for gene if not present
        exp_dict[tf][gene].append([row['time'],
                                   row['log2_cleaned_ratio']])
    return exp_dict

def sort_genes(path, sep):
    '''
    Returns the list of TFs and targets to be included in GRN based on information found in an input file. 
    Also populates thresh_df, a global dictionary containing all time and amplitude thresholds for each gene.

    Parameters:
    path (string): path of the input file
    sep (string): separator to parse input file

    Returns:
    tfs (list of strings): list of TFs
    targets (list of strings): list of target genes
    '''
    df = pd.read_csv(path, sep=sep)
    tfs = []
    targets = []
    for _, row in df.iterrows():
        gene = row['gene']
        t_thresh = row['t_thresh']
        amp_thresh = row['amp_thresh']
        thresh_df.loc[len(thresh_df.index)] = [gene, t_thresh, amp_thresh]
        gene_type = row['type']
        if 'target' in gene_type:
            targets.append(gene)
        if 'tf' in gene_type:
            tfs.append(gene)
    return tfs, targets

def filter_df(path, gene_path=None):
    '''
    Read input files to store data of desired TFs and targets into a pandas DataFrame.

    Parameters: 
    path (string): input file containing all gene expression data of all TFs and targets
    gene_path (string): input file containing TFs and targets

    Returns: 
    df (pandas DataFrame): contains gene expression data of TFs and targets
    '''
    # reading data
    df = pd.read_csv(path, sep='\t')

    # restricting nodes
    if gene_path != None:
        tfs, targets = sort_genes(gene_path, ',')
        df = df[df['TF'].isin(tfs) & df['GeneName'].isin(targets)]
        if df.empty:
            raise Exception("DataFrame is empty. Please check the input TFs and target genes.")
    
    # filter columns
    df = df[['TF', 'GeneName', 'time', 'log2_cleaned_ratio']]
    return df

def read_peak_times():
    '''
    Read peak expression times of genes to plot nodes on a timeline.

    Returns:
    peak_times (pandas DataFrame): peak times of present genes
    '''
    df = pd.read_csv('peak_times.csv')
    peak_times = df[df['Genes'].isin(gene_nodes)]
    return peak_times

# using sigmoidal fit to retrieve gene expression info
def get_sig_info(tf, gene):
    '''
    Retrieve whether the interaction is activation/inhibition and when it occurs based on sigmoidal curve fit. 
    This is based on the t1/2 of the sigmoidal curve.
    
    Parameters:
    tf (string): name of TF
    gene (string): name of target gene

    Returns:
    time of activation/inhibition (float)
    whether its activation or not (boolean)
    '''
    params = sigmoid_curve_fit(tf, gene) # [L_opt, x0_opt, k_opt, b_opt] or [L1, x01, k1, b1, L2, x02, k2, b2]
    amp_thresh = thresh_df.loc[thresh_df['gene'] == tf, 'amp_thresh'].values[0]
    if params is None or len(params) == 0 or abs(params[0]) < amp_thresh:
        return None, None
    return params[1], (params[0] > 0)

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return y

def double_sigmoid(x, L1, x01, k1, b1, L2, x02, k2, b2):
    y = (L1 / (1 + np.exp(-k1 * (x - x01))) + b1) + (L2 / (1 + np.exp(-k2 * (x - x02))) + b2)
    return y

def scale_data(data):
    '''
    Scale the data to a range suitable for clustering.

    Parameters: 
    data (list of floats): expression level data for one TF-target pair

    Returns:
    (list of floats): rescaled version of expression level data
    '''
    data_min = min(data)
    data_max = max(data)
    if data_max == data_min:
        return [0 for _ in data]
    return [(d - data_min) / (data_max - data_min) for d in data]


def sigmoid_curve_fit(tf, gene):
    '''
    Find the sigmoidal curve of best fit (single or double).

    Parameters:
    tf (string): name of TF
    gene (string): name of target gene

    Returns: 
    params (list): sigmoidal curve constants [L ,x0, k, b]
    '''
    time_series = exp_dict[tf][gene]
    xdata = []
    ydata = []
    for time in time_series:
        xdata.append(time[0])
        ydata.append(time[1])
    p0 = [max(ydata), np.median(xdata), 1, min(ydata)] 

    # prematurely remove time series data that doesn't have a big enough amplitude
    amp_thresh = thresh_df.loc[thresh_df['gene'] == tf, 'amp_thresh'].values[0]
    if max(ydata) - min(ydata) < amp_thresh:
        return None 
    
    # curve fit
    try:
        bounds = ([-np.inf, 0, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf]) # bounds and method specified
        params, _ = scipy.optimize.curve_fit(sigmoid, xdata, ydata, p0=p0, bounds=bounds, method='dogbox', maxfev=10000)
        return params
    except Exception:
        pass
    try:
        bounds = ([-np.inf, 0, 0, -np.inf, -np.inf, 0, 0, -np.inf], 
                  [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])
        params, _ = scipy.optimize.curve_fit(double_sigmoid, xdata, ydata, p0=p0, bounds=bounds, method='dogbox', maxfev=10000)
        return params
    except Exception as e:
        print(f"{e} - TF: {tf}, target: {gene}")
        return None

# building heat maps
def heat_map_colors():
    '''
    Set the color scheme of heat map to Haase lab colors.
    
    Returns:
    haase (matplotlib.colors.LinearSegmentedColormap): color scheme
    '''
    norm = matplotlib.colors.Normalize(-1.5,1.5)
    colors = [[norm(-1.5), "cyan"],
          [norm(0), "black"],
         [norm(1.5), "yellow"]]
    haase = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    return haase

def create_heat_maps(df):
    '''
    Build heat maps based on gene expression levels, and creates .png file containing output. 
    Arrange subplots based on TF. 

    Parameters:
    df (pandas DataFrame): gene expression data of desired TFs and targets
    '''
    # run heat_map_colors() to get haase color scheme
    haase = heat_map_colors()

    # num_cols and num_rows to determine dimensions of plot and subplots
    tfs = df['TF'].unique()
    num_tfs = len(tfs)
    num_cols = 4
    num_rows = math.ceil(num_tfs / num_cols)

    # create plot
    fig = plt.figure(figsize = (15,10 + num_rows * 3))
    fig.subplots_adjust(hspace=0.4, wspace=0.4, top = 0.90)
    fig.suptitle("IDEA Dataset Expression Levels", fontsize = 15)
    

    for i, tf in enumerate(tfs):    
        # Filter data for the specific TF
        tf_data = df[df['TF'] == tf]

        # Create a pivot table for heatmap 
        heatmap_data = tf_data.pivot(index='GeneName', columns='time', values='log2_cleaned_ratio')

        # Make subplot
        ax = plt.subplot(num_rows, num_cols, i + 1)
        sns.heatmap(heatmap_data, cmap=haase, cbar=True, vmin=-2, vmax=2, ax=ax)
        ax.set_title(tf)
        ax.set_xlabel('time (min)')
        ax.set_ylabel('')  

        # set tick positions and labels
        # times = np.array([0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 90.0])
        gene_names = heatmap_data.index.to_numpy()
        # ax.set_xticks(np.arange(len(times)), labels=times)
        ax.set_yticks(np.arange(len(gene_names)) + 0.25, labels=gene_names, fontsize=6)
        plt.xticks(rotation=45)  

    plt.savefig("heat_maps.png")
    plt.show()

# grn network construction
def build_ref_network():
    '''
    Build cell cycle reference network on top of generated GRN. 
    '''
    df = pd.read_csv('ref_edges.csv')
    for _, row in df.iterrows():
        tf = row['reg']
        if tf not in gene_nodes:
            gene_nodes[tf] = structs.GeneNode(tf)
        tf_node = gene_nodes[tf]

        target = row['target']
        if target not in gene_nodes:
            gene_nodes[target] = structs.GeneNode(target)
        target_node = gene_nodes[target]

        act = row['type'] == 'act'

        edge = tf_node.get_edge(target_node, act)
        if edge is None:
            tf_node.add_edge(target_node, act)


def build_tree(tf, edges_df):
    '''
    Create a tree structure with a height of 1 for a given TF.

    Parameters:
    tf (string): TF
    edges_df (pandas DataFrame): running DataFrame of current Edges
    '''
    # get time threshold
    t_thresh = thresh_df.loc[thresh_df['gene'] == tf, 't_thresh'].values[0]

    # initializing the root node if it doesn't exist already
    if tf not in gene_nodes:
        gene_nodes[tf] = structs.GeneNode(tf)
    node = gene_nodes[tf]

    # establish edges  
    for gene in exp_dict[tf].keys(): # target
        time, sign = get_sig_info(tf, gene)
        if time is not None and time > 0 and time < t_thresh:
            if gene not in gene_nodes:
                gene_nodes[gene] = structs.GeneNode(gene)
            target = gene_nodes[gene]
            node.add_edge(target, sign, time) 
            edges_df.loc[len(edges_df.index)] = [tf, gene, 'act' if sign else 'rep', time]
        
def visualize_gene_network():
    '''
    Create visualization of GRN using graphviz package. 
    Output stored in a .pdf file.
    '''
    # produce ref graph
    build_ref_network()

    # create digraph
    dot = graphviz.Digraph(comment='Gene Regulatory Network')
    
    # Add nodes
    for gene in gene_nodes:
        dot.node(gene_nodes[gene].gene, shape='box')
    
    # Add edges
    both, refs, grns = network_edges(dot)

    # create legend
    create_legend(dot)

    edge_counts_text = (
        f"Number of common edges: {len(both)}\n"
        f"Number of edges in only ref: {len(refs)}\n"
        f"Number of edges in only grn: {len(grns)}"
    )
    dot.node('edge_counts', label=edge_counts_text, shape='plaintext', fontsize='12')
    # dot.edge('edge_counts', 'main_annotation', style='invis')  # Link text to main annotation invisibly

    # Render the graph
    dot.render(f'grns/gene_network_{timestamp}', view=True)

def create_legend(dot):
    '''
    Create a legend on a given Digraph in graphviz.

    Parameters: 
    dot (Digraph): the given Digraph in graphviz
    '''
    with dot.subgraph(name='cluster_legend') as legend:
        legend.attr(label='Legend', labelloc='t', fontsize='20')
    
        legend.node('both_legend', 'In both graphs', shape='plaintext')
        legend.node('refs_legend', 'Only found in reference graph', shape='plaintext')
        legend.node('grns_legend', 'Only found in algorithmic graph', shape='plaintext')
        
        # Creating colored edges for legend
        legend.node('legend_space1', '', width='0.1', shape='plaintext')
        legend.node('legend_space2', '', width='0.1', shape='plaintext')
        legend.node('legend_space3', '', width='0.1', shape='plaintext')
        
        legend.edge('both_legend', 'legend_space1', color='#008000', arrowhead='none')
        legend.edge('refs_legend', 'legend_space2', color='#0059b3', arrowhead='none')
        legend.edge('grns_legend', 'legend_space3', color='#b30000', arrowhead='none')

def network_edges(dot):
    '''
    Create edges in a given Digraph.

    Parameters:
    dot (Digraph): the given Digraph in graphviz
    '''
    both, refs, grns = compare_edges('ref_edges.csv', f'grn_edges/grn_edges_{timestamp}.csv')
    for gene in gene_nodes:
        for edge in gene_nodes[gene].edges:
            # determine arrowhead
            if edge.act:
                arrowhead = 'normal'
            else:
                arrowhead = 'tee'
            # determine edge color
            if edge in both:
                color = '#008000'
            elif edge in refs:
                color = '#0059b3'
            elif edge in grns:
                color = '#b30000'
            else:
                color = 'black'
            dot.edge(gene, edge.target.gene, color=color, arrowhead=arrowhead)
    return both, refs, grns

def build_network(df):
    '''
    The primary function for building the overall network.

    Parameters:
    df (pandas DataFrame): contains gene expression data of desired TFs and targets
    '''
    # build network
    tfs = df['TF'].unique()
    edges_df = pd.DataFrame(columns=['reg', 'target', 'type', 'time'])
    for tf in tfs:
        build_tree(tf, edges_df) 
    edges_df.to_csv(f"grn_edges/grn_edges_{timestamp}.csv", index=False)

    visualize_gene_network()

# find paths
def find_paths(tf, target, graph, path_limit=3):
    '''
    Find all paths from one node to another.

    Parameters:
    tf (string): name of starting node
    target (string): name of ending node
    graph (string): the graph to get edges from
    path_limit (int): maximum length of paths found

    Returns: 
    paths (list of lists of Edges): each path is represented by a list of Edges. 
        All finished paths must have the target node as the target of the last Edge.
    '''
    tf_node = gene_nodes[tf]
    paths = []
    edges = get_edges_info(tf_node.gene, graph)
    for edge in tf_node.edges:
        str_edge = [tf, edge.target.gene]
        if str_edge in edges:
            new_paths = find_paths_recursive([[edge]], tf, target, graph, 1, path_limit)
            paths.extend(new_paths)
    if len(paths) == 0:
        return None
    return paths

def find_paths_recursive(paths, tf, target, graph, path_length, path_limit):
    '''
    Recursive helper function for find_paths. Tracks path recursively.

    Parameters:
    paths (list of lists of Edges): current paths that are being tracked. 
        Contains paths in progress and finished paths.
    target (string): string of the target node
    graph (string): graph to get edges from. 'ref' if reference, 'grn' if from generated GRN.
    path_length (int): current length of each path in progress in paths (they should all be the same length)
    path_limit (int): maximum length of a path
    '''
    # edge case. theoretically this shouldn't happen
    if path_length > path_limit:
        return
    
    # finish up the recursion
    if path_length == path_limit:
        final_paths = []
        for path in paths:
            last_edge = path[len(path)-1]
            last_node = last_edge.target.gene
            if last_node == target:
                final_paths.append(path)
        return final_paths
    
    # recursion
    new_paths = []
    for path in paths:
        last_edge = path[len(path)-1]
        last_node = last_edge.target # GeneNode
        last_node_name = last_node.gene
        if last_node_name == target:
            new_paths.append(path)
        elif len(last_node.edges) == 0:
            continue
        else:
            reg = last_node.gene
            edges = get_edges_info(reg, graph)
            for edge in last_node.edges:
                tar = edge.target.gene # string
                node_names = path_to_strings(tf, path)
                if tar in node_names:
                    break
                row = [reg, tar]
                if row in edges:
                    new_path = path + [edge]
                    new_paths.append(new_path)
    return find_paths_recursive(new_paths, tf, target, graph, path_length + 1, path_limit)

def get_edges_info(reg, graph):
    '''
    Given a regulator and the specified graph, find all edges from that regulator in the graph.

    Parameters:
    reg (string): name of the regulator
    graph (string): 'ref' if reference graph, 'grn' if GRN graph

    Returns:
    rows (list of lists of strings): the rows from edges_info_{timestamp}.csv that satisfy the attributes
    '''
    df = pd.read_csv(f"edges_info/edges_info_{timestamp}.csv")
    df = df[(df['reg'] == reg) & ((df['group'] == 'both') | (df['group'] == f'{graph}s'))]
    rows = []
    for _, row in df.iterrows():
        new_row = [row['reg'], row['target']]
        rows.append(new_row)
    return rows


def paths_to_df(graph, path_limit=3):
    '''
    Outputs all paths into a DataFrame to be converted to a .csv file.

    Parameters:
    graph (string): specifies which graph to get paths from. 
        'ref' if reference graph, 'grn' if generated GRN.
    path_limit (int): the maximum length of every path

    Returns: 
    path_df (pandas DataFrame): DataFrame to be 
    '''
    # getting proper edges
    if graph != 'ref' and graph != 'grn':
        raise ValueError("Improper input for graph attribute.")
    
    # getting paths
    path_df = pd.DataFrame(columns=['start', 'end', 'path', 'graph', 'sign', 'time'])
    genes = list(gene_nodes.keys())
    for i in range(len(genes)):
        for j in range(len(genes)):
            gene1 = genes[i]
            gene2 = genes[j]
            paths = find_paths(gene1, gene2, graph, path_limit)
            if paths is None:
                continue
            else:
                for p in paths:
                    path = path_to_strings(gene1, p)
                    if graph == 'grn':
                        time = path_time(p)
                    else:
                        time = np.nan
                    sign = path_sign(p)
                    path_df.loc[len(path_df.index)] = [gene1, gene2, path, graph, sign, time]
    print(path_df)
    return path_df

def path_time(path):
    '''
    Calculate the total amount of time it takes for a full path to activate/repress.

    Parameters:
    path (list of Edges): represents the Edges in a path

    Returns:
    time (float): total time of path
    '''
    time = 0.0
    for edge in path:
        t = edge.time
        if t is None:
            return np.nan
        time = time + t
    return time

def path_sign(path):
    '''
    Calculate the sign (activation/repression) of a given path.

    Parameters:
    path (list of Edges): represents the Edges in a path

    Returns:
    sign (string): 'act' if activating, 'rep' if repressing
    '''
    if path is None or len(path) == 0:
        return None
    neg_edges = 0 # negative edge counter
    for edge in path:
        if not edge.act:
            neg_edges += 1
    if neg_edges % 2 == 1:
        return 'rep'
    return 'act'

def find_all_paths(path_limit=3):
    '''
    Convert path DataFrames into a output .csv file.
    '''
    df1 = paths_to_df('ref', path_limit=path_limit)
    df2 = paths_to_df('grn', path_limit=path_limit)
    print(len(df1.index))
    print(len(df2.index))
    df = pd.concat([df1, df2], axis=0).reset_index(drop=True)
    df = df.sort_values(by=['graph', 'start', 'end']).reset_index(drop=True)
    df.to_csv("paths.csv")

def path_to_strings(tf, path):
    '''
    Given a path, return a list of the names of the nodes that the path traverses through.

    Parameters:
    tf (string): string representing the first node
    path (list of Edges): path to be converted to list of strings

    Returns:
    strings (list of strings): list of names of nodes in path
    '''
    strings = []
    if tf is not None:
        strings = [tf]
    for edge in path:
        target = edge.target
        node = target.gene
        strings.append(node)
    return strings

# compare with reference graph
def compare_edges(ref, grn):
    '''
    Compare edges in ref_edges.csv and grn_edges.csv to find agreements and discrepancies. 
    Outputs a .csv file of edges (i.e. reg and target) and which list they belong in.

    Parameters:
    ref (string): path to ref_edges.csv
    grn (string): path to grn_edges.csv

    Returns: 
    both (list of Edges): Edges found in both files
    refs (list of Edges): Edges found only in ref_edges.csv and not in grn_edges.csv
    grns (list of Edges): Edges found only in grn_edges.csv and not in ref_edges.csv
    '''
    # read csv files into DataFrames
    ref_df = pd.read_csv(ref)
    grn_df = pd.read_csv(grn)
    grn_df = grn_df.loc[:, grn_df.columns.difference(['time'])]

    ref_edges = set(tuple(row) for row in ref_df.to_records(index=False))
    grn_edges = set(tuple(row) for row in grn_df.to_records(index=False))

    # Find common edges
    both = ref_edges & grn_edges

    # Find edges only in ref_edges.csv
    refs = ref_edges - grn_edges

    # Find edges only in grn_edges.csv
    grns = grn_edges - ref_edges

    # Convert sets back to lists of tuples
    both = list(both)
    refs = list(refs)
    grns = list(grns)

    # DataFrame of all edges to be converted to .csv file
    all_edges = pd.DataFrame(columns=['reg', 'target', 'act', 'group'])

    # convert to Edges and add to edges_df
    for i in range(len(both)):
        temp = both[i]
        reg = temp[0]
        target = temp[1]
        act = temp[2] == 'act'
        both[i] = gene_nodes[reg].get_edge(gene_nodes[target], act)
        all_edges.loc[len(all_edges.index)] = [reg, target, act, 'both']

    for i in range(len(refs)):
        temp = refs[i]
        reg = temp[0]
        target = temp[1]
        act = temp[2] == 'act'
        refs[i] = gene_nodes[reg].get_edge(gene_nodes[target], act)
        all_edges.loc[len(all_edges.index)] = [reg, target, act, 'refs']
    
    for i in range(len(grns)):
        temp = grns[i]
        reg = temp[0]
        target = temp[1]
        act = temp[2] == 'act'
        grns[i] = gene_nodes[reg].get_edge(gene_nodes[target], act)
        all_edges.loc[len(all_edges.index)] = [reg, target, act, 'grns']

    # render .csv file
    all_edges.to_csv(f"edges_info/edges_info_{timestamp}.csv")

    # return lists of Edges
    return both, refs, grns

def df_to_edges(df):
    '''
    Given a pandas DataFrame containing edge information, return the Edges listed.

    Parameters: 
    df (pandas DataFrame): input DataFrame with columns 'reg', 'target', and 'type'

    Returns:
    edges (list of Edges): the list of Edges corresponding to the input df
    '''
    edges = []
    for _, row in df.iterrows():
        reg = gene_nodes[row['reg']] if row['reg'] in gene_nodes else None
        target = gene_nodes[row['target']] if row['target'] in gene_nodes else None
        if reg == None or target == None:
            break
        act = row['type'] == 'act'
        edge = reg.get_edge(target, act)
        edges.append(edge)
    return edges

# main function
def main():
    path = os.path.join(DATA_DIR, DATA_FILE)
    gene_path = os.path.join(GENE_DIR, GENE_FILE)

    df = filter_df(path, gene_path)

    extract_expression_data(df)

    thresh_df = pd.read_csv(gene_path)
    thresh_df = thresh_df.columns.difference(['type'])

    # # create heat maps
    # create_heat_maps(df)

    # reformat data and build network
    build_network(df)
    find_all_paths()
    
if __name__ == '__main__':
    main()   