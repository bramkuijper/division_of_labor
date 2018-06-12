#!/usr/bin/env python3

# plot the division of labor simulations

import pandas as pd
import itertools
import subprocess 
import argparse
import numpy as np
import sys
import matplotlib
matplotlib.use("Agg")

from pathlib import Path
from pathlib import PurePath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm

plt.style.use('base')

# some stuff to render fonts in graphs
rcParams['axes.labelsize'] = 15
rcParams['text.usetex'] = True
rcParams['font.family'] = 'sans-serif'

# some stuff to render fonts in graphs
# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]


# basename of the executable that makes histgrams out of the allele frequency
# data. Such histograms are necessary for plotting things
histo_exe_basename = "xreadhisto"

# the basename of the resulting histogram file
histogram_basename = "histograms.csv"
number_columns_histo = 6

# filename of the work allocation file
work_file_basename = "data_work_alloc.txt"

# the file with the allelic distributions (we need to make a histogram
# of this) 
allele_file_basename = "threshold_distribution.txt"

#########################################

#           read in the data

#########################################

# first set up argument parser object
# so that we can pass some command line flags to the python script
parser = argparse.ArgumentParser(description="Plot simulation result for the fixed response threhsold model.")

# add the data_work_alloc file as an obligatory argument
parser.add_argument('pathname', 
        metavar="directory",
        type=str,
        nargs=1,
        help="a pathname in which the " + work_file_basename + " file is located")

parser.add_argument('--hist', 
        dest='use_histogram', 
        action='store_const',
        const=True,
        default=False,
        help="specify a histogram file to plot as well")

# ok parse the command line arguments
args = parser.parse_args()

# resolve the current path (e.g., prevent symlinks etc)
current_simpath = Path(vars(args)["pathname"][0]).resolve()

# now get the file name of the work alloc file within this path
work_alloc_file = list(current_simpath.glob(work_file_basename))

assert(len(work_alloc_file) == 1)

# also store the absolute base folder in which the data_work_alloc.txt is contained
base_folder = current_simpath

# whether we need to produe the histograms showing branching in the 
# different trait values
use_hist = vars(args)["use_histogram"]

# get the work allocation data and plot it
work_alloc_data = pd.read_csv(work_alloc_file[0],
        sep="\t")

# ok we need to plot histograms of the allelic values
# hence we need to process the allelic values data
if use_hist:
    # perform a system call to generate histogram from the allele frequency data
    # https://stackoverflow.com/questions/3718657/how-to-properly-determine-current-script-directory 
    script_dir = PurePath(Path(__file__).resolve())

    # generate file name of the histogram generator executable
    histo_exe = Path(script_dir.parent / histo_exe_basename)

    # see whether the exe is indeed there
    assert(histo_exe.exists())

    # also check if the histo name is already existing, in which case
    # we do not have to run the thingy again
    result_histo_name = base_folder / histogram_basename

    if (not result_histo_name.exists()):

        # then create the filename of the allele_distrib_1.txt containing
        # all the allelic values
        allele_dist_file = base_folder / allele_file_basename

        assert(allele_dist_file.exists())

        # call the histogram creation executable and let it work on the file
        subprocess.call([str(histo_exe), 
            str(allele_dist_file), 
            str(number_columns_histo)])

        # see whether the histograms.csv file is produced
        assert(result_histo_name.exists())

    # get data on branching for learning and forgetting
    # this is a histogram with all sorts of values
    branch_data = pd.read_csv(
            filepath_or_buffer=result_histo_name, 
            sep=";", 
            index_col=False
            )

    print("read branch data")

    # now 4th root transform the data, so that also see rare frequencies
    def the_pow(x):
        return(x**(0.25))

    # power transform the count data
    branch_data["count"] = branch_data[
            "count"].apply(the_pow)

    # generate the pivot table that is necessary to plot it in imshow()
    def generate_pivot(the_data, x, y, z):

        # make a pivot table
        the_pivot = the_data.pivot_table(
                values=z, 
                index=y, 
                columns=x)

        x, y = np.meshgrid(
                the_pivot.columns.values, 
                the_pivot.index.values)

        z = the_pivot.values

        return(x, y, z)


    # generate a pivot table for the female threshold data
    (x_threshold1, y_threshold1, threshold1_count) = generate_pivot(
            the_data = branch_data[branch_data["traitname"]=="trait0"], 
            x="generation",
            y="bin_start",
            z="count"
            )

    print("pivot 1")

    # generate a pivot table for the female threshold data
    (x_threshold2, y_threshold2, threshold2_count) = generate_pivot(
            the_data = branch_data[branch_data["traitname"]=="trait1"], 
            x="generation",
            y="bin_start",
            z="count"
            )

#########################################

#           make figure

#########################################

# initialize the figure
fig = plt.figure(figsize=(10,18))

# generate the grid of the graph
# see: 
widths = [ 1, 0.05 ]

numrows = 5

heights = [ 1 for x in range(0,numrows) ]

# make the grid
gs = gridspec.GridSpec(
        nrows=len(heights),
        ncols=len(widths),
        width_ratios=widths,
        height_ratios=heights)

if use_hist:
    ax = plt.subplot(gs[0,0])

    # the plot for learning
    ax.imshow(threshold1_count,
        cmap="jet",
        extent=[x_threshold1.min(), 
            x_threshold1.max(), 
            y_threshold1.min(), 
            y_threshold1.max()],
        origin="lower",
        aspect="auto")

    ax.set_ylabel(r"Threshold task 1")

    # start next entry of the graph
    ax = plt.subplot(gs[1,0])

    # the plot for learning
    ax.imshow(threshold2_count,
        cmap="jet",
        extent=[x_threshold2.min(), 
            x_threshold2.max(), 
            y_threshold2.min(), 
            y_threshold2.max()],
        origin="lower",
        aspect="auto")

    ax.set_ylabel(r"Threshold task 2")


# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[2,0])

# ok, now calculate min, max, mean, stddev and median
# for all variables
work_alloc_data_agg = work_alloc_data.groupby("Gen").agg(['min','max','mean','std','median'])


# see https://stackoverflow.com/questions/41809118/get-columns-from-multiindex-dataframe-with-named-labels 
work_alloc_data_agg = work_alloc_data_agg.sort_index(axis=1)
idx = pd.IndexSlice
print(len(list(work_alloc_data_agg.index.get_level_values("Gen"))))
print(work_alloc_data.head())
print(work_alloc_data_agg.head())

# plot work allocation
ax = plt.subplot(gs[2,0])


# auxiliary function to get confidence envelope
def confidence_interval(variable_name):
    global work_alloc_data_agg
    
    std_min = work_alloc_data_agg[variable_name]["mean"] - work_alloc_data_agg[variable_name]["std"]
    std_max = work_alloc_data_agg[variable_name]["mean"] + work_alloc_data_agg[variable_name]["std"]

    return(std_min,std_max)



(min_sd1, max_sd1) = confidence_interval("WorkAlloc1")
(min_sd2, max_sd2) = confidence_interval("WorkAlloc1")

# calculate lower bound of confidence envelope by subtracting standard
# deviation from the mean


ax.plot(
        work_alloc_data_agg.loc[:, idx["WorkAlloc1",["mean"]]],
        color="blue",
        label="Work alloc task 1")

ax.plot(
        min_sd1,
        color="lightblue",
        label="_nolabel")

ax.plot(
        max_sd1,
        color="lightblue",
        label="_nolabel")


ax.plot(
        work_alloc_data_agg.loc[:, idx["WorkAlloc2",["mean"]]],
        color="red",
        label="Work alloc task 2")


ax.plot(
        min_sd2,
        color="pink",
        label="_nolabel")

ax.plot(
        max_sd2,
        color="pink",
        label="_nolabel")

# add a legend
ax.legend()

ax.set_ylabel(r"Work alloc")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)


(min_sd_idle, max_sd_idle) = confidence_interval("Idle")

# start next entry of the graph
ax = plt.subplot(gs[3,0])

ax.plot(
        work_alloc_data_agg.loc[:, idx["Idle",["mean"]]],
        label="Idle",
        color="blue")

ax.plot(
        min_sd_idle,
        color="lightblue",
        label="_nolabel")

ax.plot(
        max_sd_idle,
        color="lightblue",
        label="_nolabel")

ax.set_ylabel(r"Idle")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)

# start next entry of the graph
ax = plt.subplot(gs[4,0])

(min_sd_fitness, max_sd_fitness) = confidence_interval("Fitness")

ax.plot(
        work_alloc_data_agg.loc[:, idx["Fitness",["mean"]]],
        color="blue")

ax.plot(
        min_sd_fitness,
        color="lightblue")

ax.plot(
        max_sd_fitness,
        color="lightblue")

ax.set_ylabel(r"Total $w$")

format = "pdf"

plt.savefig("graph_simplot." + format, 
        format=format, 
        bbox_inches="tight")
