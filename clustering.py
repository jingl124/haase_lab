import os
import datetime
import seaborn as sns
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import grn_finder as grn
import math

# Global timestamp
timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# Constants
DATA_DIR = "/Users/jingliu/Documents/haase/IDEA_data"
DATA_FILE = "idea_tall_expression_data.tsv"
PLOT_OUTPUT_DIR = f'cluster_plots/cluster_plots_{timestamp}'
HEAT_OUTPUT_DIR = f'heat_map_clusters/heat_maps_{timestamp}'

# Utility Functions
def scale_data(data):
    """
    Scale the data to a range suitable for clustering.

    Parameters: 
    data (list of floats): expression level data for one TF-target pair

    Returns:
    (list of floats): rescaled version of expression level data
    """
    data_min = min(data)
    data_max = max(data)
    if data_max == data_min:
        return [0 for _ in data]
    return [(d - data_min) / (data_max - data_min) for d in data]

def extract_features(series):
    """
    Extracts features from a time series.

    Parameters:
    series (pd.Series): Time series data

    Returns:
    dict: Extracted features
    """
    peaks, _ = scipy.signal.find_peaks(series)
    valleys, _ = scipy.signal.find_peaks(series * -1)
    slope = np.gradient(series)
    curvature = np.gradient(slope)
    scaled_series_np = series.to_numpy() if isinstance(series, pd.Series) else np.array(series)
    fft_values = np.abs(scipy.fft.fft(scaled_series_np))
    dominant_freq = np.argmax(fft_values[1:len(fft_values)//2]) + 1  # Ignore the zero frequency

    features = {
        'num_peaks': len(peaks),
        'num_valleys': len(valleys),
        'mean_slope': np.mean(slope),
        'mean_curvature': np.mean(curvature),
        'dominant_freq': dominant_freq,
        'skewness': scipy.stats.skew(series),
    }
    return features

def scale_features():
    """
    Scales the features of the time series data.

    Returns:
    tuple: DataFrame of features, scaled features array
    """
    feature_list = []
    for tf in grn.exp_dict:
        for target in grn.exp_dict[tf]:
            time_series_data = grn.exp_dict[tf][target]
            times, values = zip(*time_series_data)
            values = scale_data(values)
            series = pd.Series(data=values, index=times)
            row = {'TF': tf, 'target': target}
            features = extract_features(series)
            row.update(features)
            feature_list.append(row)
    feature_df = pd.DataFrame(feature_list)

    # Scale features
    columns_to_scale = feature_df.columns.difference(['TF', 'target'])
    scaler = StandardScaler()
    scaled_df = feature_df[columns_to_scale]
    features_scaled = scaler.fit_transform(scaled_df)

    return feature_df, features_scaled

# Clustering Functions
def elbow_curve():
    """
    Plots the elbow curve to determine the optimal number of clusters.

    Parameters:
    features_scaled (array): Scaled features
    """
    _, features_scaled = scale_features()
    wcss = []
    max_clusters = 20
    for k in range(1, max_clusters + 1):
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(features_scaled)
        wcss.append(kmeans.inertia_)

    plt.figure(figsize=(8, 5))
    plt.plot(range(1, max_clusters + 1), wcss, marker='o')
    plt.title('Elbow Method for Optimal k')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Within-Cluster Sum of Squares (WCSS)')
    plt.savefig("elbow.png")

def classify_time_series(num_clusters):
    """
    Applies KMeans clustering to classify time series data.

    Parameters:
    num_clusters (int): Number of clusters

    Returns:
    DataFrame: DataFrame with cluster labels
    """
    feature_df, features_scaled = scale_features()
    kmeans = KMeans(n_clusters=num_clusters, random_state=42)
    clusters = kmeans.fit_predict(features_scaled)
    feature_df['cluster'] = clusters
    columns = ['TF', 'target', 'cluster']
    cluster_df = feature_df[columns]
    
    return cluster_df

def plot_clusters(feature_df, scaled=False):
    """
    Plots time series for each cluster.

    Parameters:
    num_clusters (int): Number of clusters
    scaled (bool): Whether to scale the data
    """
    # feature_df = classify_time_series(num_clusters)
    # feature_df.to_csv("clusters.csv")
    # os.makedirs(OUTPUT_DIR, exist_ok=True)
    num_clusters = feature_df['cluster'].max() + 1

    for cluster in range(num_clusters):
        cluster_df = feature_df[feature_df['cluster'] == cluster]
        plt.figure(figsize=(10, 6))
        for _, row in cluster_df.head().iterrows():  # Plot first 5 series for brevity
            tf = row['TF']
            target = row['target']
            time_series = grn.exp_dict[tf][target]
            timestamps, values = zip(*time_series)
            if not scaled:
                plt.plot(timestamps, values, marker='o', linestyle='-', label=f'{tf}-{target}')
            else: 
                plt.plot(timestamps, scale_data(values), marker='o', linestyle='-', label=f'{tf}-{target}')
        plt.title(f'Cluster {cluster}')
        plt.legend()
        os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)
        plt.savefig(os.path.join(PLOT_OUTPUT_DIR, f'cluster_{cluster}.png'))
        plt.close()

def cluster_heat_maps(df, cluster_df):
    '''
    Plot out gene expression heat maps by cluster.
    '''
    haase = grn.heat_map_colors()

    # num_cols and num_rows to determine dimensions of plot and subplots
    num_clusters = cluster_df['cluster'].max() + 1
    num_cols = 4
    num_rows = math.ceil(num_clusters / num_cols)

    # create plot
    fig = plt.figure(figsize = (15,10 + num_rows * 3))
    fig.subplots_adjust(hspace=0.4, wspace=0.4, top = 0.90)
    fig.suptitle("Expression Levels by Cluster", fontsize = 15)

    for i in range(num_clusters):
        temp_df = cluster_df[cluster_df['cluster'] == i]
        for _, row in temp_df.iterrows():
            tf = row['TF']
            target = row['target']           
            data = df[df['TF'] == tf & df['target'] == target]

            # Create a pivot table for heatmap 
            heatmap_data = data.pivot(index='GeneName', columns='time', values='log2_cleaned_ratio')

            if combined_heatmap_data.empty:
                combined_heatmap_data = heatmap_data
            else:
                combined_heatmap_data = pd.concat([combined_heatmap_data, heatmap_data])

        # Make subplot
        ax = plt.subplot(num_rows, num_cols, i + 1)
        sns.heatmap(combined_heatmap_data, cmap=haase, cbar=True, vmin=-2, vmax=2, ax=ax)
        ax.set_title(f'Cluster {i}')
        ax.set_xlabel('time (min)')
        ax.set_ylabel('')
    os.makedirs(HEAT_OUTPUT_DIR, exist_ok=True)
    plt.savefig(os.path.join(HEAT_OUTPUT_DIR, f"heatmap_clusters_{timestamp}.png"))

# Main Function
def main():
    print("Running clustering")
    path = os.path.join(DATA_DIR, DATA_FILE)
    df = grn.filter_df(path)
    grn.extract_expression_data(df)
    # elbow_curve()
    cluster_df = classify_time_series(18)
    cluster_df.to_csv("clusters.csv")
    plot_clusters(cluster_df)
    cluster_heat_maps(df, cluster_df)
    print("finished running")

if __name__ == '__main__':
    main()
