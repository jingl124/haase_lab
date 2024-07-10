# clustering the entire IDEA dataset (time series curves)
import pandas as pd
import os
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import datetime
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import antropy as ant  # For entropy calculation
import statsmodels.tsa.stattools # import acf, pacf
import grn_finder as grn

# global timestamp
timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# clustering
def extract_features(series):
    peaks, _ = scipy.signal.find_peaks(series)
    valleys, _ = scipy.signal.find_peaks(series * -1)
    slope = np.gradient(series)
    curvature = np.gradient(slope)
    scaled_series_np = series.to_numpy() if isinstance(series, pd.Series) else np.array(series)
    fft_values = np.abs(scipy.fft.fft(scaled_series_np))
    dominant_freq = np.argmax(fft_values[1:len(fft_values)//2]) + 1  # Ignore the zero frequency

    # Autocorrelation and partial autocorrelation
    autocorr = statsmodels.tsa.stattools.acf(series, nlags=10)
    # partial_autocorr = statsmodels.tsa.stattools.pacf(series, nlags=10)
    
    # Entropy
    entropy = ant.perm_entropy(series, normalize=True)
    
    # Hurst Exponent
    hurst_exponent = ant.hjorth_params(series)[0]  # Hjorth mobility can be used as a proxy for Hurst exponent
    
    # Energy
    energy = np.sum(series**2)
    
    # Root Mean Square (RMS)
    rms = np.sqrt(np.mean(series**2))
    
    # Zero-Crossing Rate
    zero_crossings = len(np.where(np.diff(np.sign(series)))[0])
    
    # Spectral Centroid and Bandwidth
    spectral_centroid = np.sum(np.arange(len(fft_values)) * fft_values) / np.sum(fft_values)
    spectral_bandwidth = np.sqrt(np.sum((np.arange(len(fft_values)) - spectral_centroid)**2 * fft_values) / np.sum(fft_values))
    
    # Area Under the Curve (AUC)
    auc = scipy.integrate.trapz(series)
    
    # Linearity (slope of linear fit)
    linear_fit = np.polyfit(np.arange(len(series)), series, 1)
    slope_of_linear_fit = linear_fit[0]

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
    print(scaled_df)
    features_scaled = scaler.fit_transform(scaled_df)

    return feature_df, features_scaled

def elbow_curve(features_scaled):
    # Elbow method to determine the optimal number of clusters
    wcss = []
    max_clusters = 10
    for k in range(1, max_clusters + 1):
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(features_scaled)
        wcss.append(kmeans.inertia_)

    # Plotting the elbow curve
    plt.figure(figsize=(8, 5))
    plt.plot(range(1, max_clusters + 1), wcss, marker='o')
    plt.title('Elbow Method for Optimal k')
    plt.xlabel('Number of Clusters (k)')
    plt.ylabel('Within-Cluster Sum of Squares (WCSS)')
    plt.show()

def classify_time_series(num_clustered):
    # Apply KMeans clustering

    kmeans = KMeans(n_clusters=num_clusters, random_state=2)
    clusters = kmeans.fit_predict(features_scaled)
    feature_df['cluster'] = clusters
    
    return feature_df

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

def plot_clusters(num_clusters, scaled=False):
    feature_df = classify_time_series(num_clusters)
    print(feature_df)
    feature_df.to_csv("clusters.csv")
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
        # plt.show()
        output_dir = f'cluster_plots/cluster_plots_{timestamp}'
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, f'cluster_{cluster}.png'))
        plt.close() 

# main function
def main():
    dir = "/Users/jingliu/Documents/haase/IDEA_data"
    file = "idea_tall_expression_data.tsv"
    path = os.path.join(dir, file)
    df = grn.filter_df(path)

    grn.extract_expression_data(df)
    
    # get clusters and plot line graph samples
    plot_clusters(8, scaled=False)
    
if __name__ == '__main__':
    main()  