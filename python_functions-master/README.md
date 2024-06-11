# Python functions for the Haase Lab

## Requirements:
Python, pandas, matplotlib, seaborn, IPython, scipy, math, numpy

## Files
**heatmap.py**: File containing heatmapping functions.<br />
**stripey.py**: File containing code for stripey detection, interpolation, and quantile normalization.

## Usage
### Heatmap:
To use heatmap.py in a jupyter notebook, run "%run " followed by the path to the file. Eg. "%run /Users/Sophie/Documents/python_functions/heatmap.py".
This will load the functions contained in the file for use. Run a function with the desired inputs to plot the heatmap.

Typical run:
1) %run /Users/Sophie/Documents/python_functions/heatmap.py
2) load data to be heatmapped
3) set up figure with desired subplots, size, title, etc
3) plot heatmap using desired function as described below

### Stripey:
To use stripey.py or interactive_heatmap.py in a jupyter notebook, run "%run " followed by the path to the file. Eg. "%run /Users/Sophie/Documents/python_functions/stripey.py".
This will load the functions contained in the file for use. Run a function with the desired inputs to plot the heatmap.

Typical run:
1) %run /Users/Sophie/Documents/python_functions/stripey.py
2) load data
3) detect stripeys, quantile normalize, and/or interpolate timepoints



## Heatmap Functions
#### heatmap.py:
**heatmap_max(data, gene_list, first_period)**

input:
1. data frame with time series data (index = gene names)
2. gene list with the genes to be plotted
3. first period (Note: can either be the integer corresponding to the timepoint name (default) or the column index (if first_period_name is set to False))

Z-score normalizes the dataframe and plots a heatmap ordered based off of the location of the maximum in the first period.
Returns the order of the genes as plotted.


**heatmap_order(data, order)**

input:
1. data frame with the time series data (index = gene names)
2. gene list in the desired order

Z-score normalizes the data frame and plots a heat map ordered based off of the ordered gene list.
Returns the order as plotted.


**heatmap_LE(data, gene_list, first_period)**

input:
1. data frame with the time series data (index = gene names)
2. gene list with the genes to be plotted
3. first period

Z-score normalizes the data frame and plots a heat map ordered based off of the left edge in the first period.
Converted partially from R function heatOrderLE_firstper.  
Returns the order of the genes as plotted.


**Optional inputs for functions:**
1. first_period_name: gives option to supply either the integer corresponding to the name of the column (when set to True) or the integer corresponding to the index of the column (when set to False)
default = True
2. yticks: gives option to include or remove y-ticks in the figure
default = False
3. cbar_bool: gives option to include or remove the color bar in the figure
default = True
4. axis: gives the option to include a subplot axis to the heatmap
default = None

#### stripey.py:
**stripey_detector(df, threshhold, df_return)**

input:
1. data frame with the time series data (index = gene names)

Optional inputs:
1. threshhold: p-value threshhold to determine cutoff for stripey determination
default = 0
2. df_return: gives option to return the stripey df with p-values and stripey calls in addition to the printed stripey list.
default = False


**interpolate_timepoints(df, timepoints, method_option)**

input:
1. data frame with the time series data (index = gene names)
2. list of timepoints to interpolate

optional input:
1. method_option: option for type of interpolation
default = "pchip"

**qn_normalize(df)**

input:
1. data frame with the time series data (index = gene names)




Author: Sophie Campione, Robert Moseley
Haase Lab, 2021
