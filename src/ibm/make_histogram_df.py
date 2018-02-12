#!/usr/bin/python3

# make a histogram out of dataframe that has
# generation;trait_value1;trait_value2;...
# in the shape of
# generation;bin_trait1;count_trait1;bin_trait2;count_trait2,...


import pandas as pd
import sys
import os
import copy
import numpy as np

nbins = 100

# get dataset
branch_data = pd.read_csv(
        filepath_or_buffer=sys.argv[1], 
        sep=";", 
        )

# get a list of column names
col_names = list(branch_data.columns.values)

# set a list of column names for the
# future histogram dataframe
col_names_hist = ["generation"]

# get a list of column names of which
# to make a histogram
yvars = copy.copy(col_names)
yvars.remove("generation")

# make empty dictionary to contain the
# minimum and maximum values for each column
# for which a histogram needs to be made
min_max_vals = {}

# loop through column names and make a '_bin'
# and a '_count' column for each column name
for col_name in col_names:

    if col_name != "generation":
        col_names_hist += [ 
                col_name + "_bin",
                col_name + "_count"
                ]

        min_val = branch_data[col_name].min()
        max_val = branch_data[col_name].max()

        min_max_vals[col_name] = (min_val,max_val)

hist_data_frame = pd.DataFrame(
        columns = col_names_hist)

# get a list of all generations
generations = list(branch_data["generation"].unique())

# now go through each generation and make histo
for generation_i in generations:

    print(generation_i)

    # initialize subdataframe which soon will be
    # merged with main dataframe
    curr_data_frame = pd.DataFrame(columns = col_names_hist)

    # fill with a lot of values of this generation
    curr_data_frame["generation"] = [
            generation_i for i in range(0,100)]

    for yvar_i in yvars:

        # make a histogram for each generation and store it
        histo = np.histogram(
                a = branch_data[yvar_i], 
                bins = nbins,
                range = min_max_vals[yvar_i])

        curr_data_frame[yvar_i + "_bin"] = pd.Series(histo[1])
        curr_data_frame[yvar_i + "_count"] = pd.Series(histo[0])

    hist_data_frame = hist_data_frame.append(
            curr_data_frame, ignore_index=True)


hist_data_frame.to_csv("hist_" + os.path.basename(sys.argv[1]),
        sep=";")
