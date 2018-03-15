#!/usr/bin/env python3

# plot the division of labor simulations

import pandas as pd
import itertools
import argparse
import numpy as np
import sys
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


# directory in which we can find the executable
# xreadhisto
histo_exe = ""


#########################################

#           read in the data

#########################################

# first set up argument parser object
# so that we can pass some command line flags to the python script
parser = argparse.ArgumentParser(description="Plot simulation result.")

parser.add_argument('--hist', 
        dest='use_histogram', 
        action='store_const',
        const=True,
        default=False)

# ok parse the command line arguments
args = parser.parse_args()

use_hist = vars(args)["use_histogram"]






# get work allocation data
# first read in the corresponding header file
f = open("header_1.txt")
fl = f.readlines()
f.close()

# the header is now one tab-separated line, e.g.,
# header1   header2     header3
# but we need to get a list of values
header = fl[0].strip().split("\t")

# get the work allocation data and plot it
work_alloc_data = pd.read_csv("data_work_alloc_1.txt",
        sep="\t",
        header=None,
        index_col=False,
        names=header) # feed the header to the csv file

# ok we need to plot histograms of the allelic values
# hence we need to process the allelic values data
if use_hist:
    # perform a system call to generate histogram from the allele frequency data
    # https://stackoverflow.com/questions/3718657/how-to-properly-determine-current-script-directory 
    script_dir = PurePath(Path(__file__).resolve())

    # generate file name of the histogram read executable
    histo_exe = Path(script_dir.parent / "xreadhisto")

    assert(histo_exe.exists())

    subprocess.call([str(histo_exe),"allele_distrib_1.txt"])


    # get data on branching for learning and forgetting
    # this is a histogram with all sorts of values
    branch_data = pd.read_csv(
            filepath_or_buffer=sys.argv[1], 
            sep=";", 
            header=None, # no header
            index_col=False,
            names=["minbin","maxbin","generation","trait","count"] # hence provide header names
            )

    print("read branch data")

    # now 4th root transform the data, so that also see rare frequencies
    def the_pow(x):
        return(x**(0.25))

    branch_data["count"] = branch_data["count"].apply(the_pow)

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


    # generate a pivot table for learning
    (x_learn, y_learn, learn_count) = generate_pivot(
            the_data = branch_data[branch_data["trait"]=="learn"], 
            x="generation",
            y="minbin",
            z="count"
            )

    print("pivot 1")

    # generate a pivot table for forgetting
    (x_forget, y_forget, forget_count) = generate_pivot(
            the_data = branch_data[branch_data["trait"]=="forget"], 
            x="generation",
            y="minbin",
            z="count"
            )

    print("pivot 2")

#########################################

#           make figure

#########################################

# initialize the figure
fig = plt.figure(figsize=(10,14))

# generate the grid of the graph
# see: 
widths = [ 1, 0.05 ]

numrows = 8

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
    ax.imshow(learn_count,
        cmap="jet",
        extent=[x_learn.min(), 
            x_learn.max(), 
            y_learn.min(), 
            y_learn.max()],
        origin="lower",
        aspect="auto")

    ax.set_ylabel(r"Learning")

    # start next entry of the graph
    ax = plt.subplot(gs[1,0])

    # the plot for learning
    ax.imshow(forget_count,
        cmap="jet",
        extent=[x_forget.min(), 
            x_forget.max(), 
            y_forget.min(), 
            y_forget.max()],
        origin="lower",
        aspect="auto")

    ax.set_ylabel(r"Forgetting")


# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[2,0])


work_alloc_data_agg = work_alloc_data.groupby("Gen").agg(['min','max','mean','std','median'])

# see https://stackoverflow.com/questions/41809118/get-columns-from-multiindex-dataframe-with-named-labels 
work_alloc_data_agg = work_alloc_data_agg.sort_index(axis=1)
idx = pd.IndexSlice
print(work_alloc_data_agg.loc[:, idx["FitWork1",["mean"]]])

work_alloc_data_agg.columns = work_alloc_data_agg.columns.droplevel(0)

print(work_alloc_data_agg.columns.describe())

sys.exit(1)


ax.plot(
        data_agg["Gen"],
        data_agg["FitWork1"],
        marker=".",
        markerfacecolor="blue",
        markeredgecolor="blue")

ax.plot(
        work_alloc_data["Gen"],
        work_alloc_data["FitWork2"],
        marker=".",
        markerfacecolor="red",
        markeredgecolor="red")

# add a legend
ax.legend(
        (r'$w_{1}$',
            r'$w_{2}$'))

ax.set_ylabel(r"Fitwork")


# start next entry of the graph
ax = plt.subplot(gs[3,0])

ax.plot(
        x=work_alloc_data["Gen"],
        y=work_alloc_data["WorkAlloc1"],
        color="blue")

ax.plot(
        x=work_alloc_data["Gen"],
        y=work_alloc_data["WorkAlloc2"],
        color="red")

# add a legend
ax.legend(
        (r'$\mathcal{W}_{1}$',
            r'$\mathcal{W}_{2}$'))


ax.set_ylabel(r"WorkAlloc")

# start next entry of the graph
ax = plt.subplot(gs[4,0])

ax.plot(
        x=work_alloc_data["Gen"],
        y=work_alloc_data["Idle"],
        color="blue")

ax.plot(
        x=work_alloc_data["Gen"],
        y=work_alloc_data["Inactive"],
        color="red")

# add a legend
ax.legend(
        (r'Idle',
            r'Inactive'))


# start next entry of the graph
ax = plt.subplot(gs[5,0])

ax.plot(
        x=work_alloc_data["Gen"],
        y=work_alloc_data["Fitness"],
        color="blue")

# start next entry of the graph
ax = plt.subplot(gs[6,0])

ax.plot(
        x=work_alloc_data["Gen"],
        y=work_alloc_data["End_stim1"],
        color="blue")

ax.plot(
        x=work_alloc_data["Gen"],
        y=work_alloc_data["End_stim2"],
        color="red")

# add a legend
ax.legend(
        (r'End_stim1',
            r'End_stim2'))


# start next entry of the graph
ax = plt.subplot(gs[7,0])

ax.plot(
        x=work_alloc_data["Gen"],
        y=work_alloc_data["mean_switches"],
        color="blue")

ax.plot(
        x=work_alloc_data["Gen"],
        y=work_alloc_data["mean_workperiods"],
        color="red")

# add a legend
ax.legend(
        (r'mean_switches',
            r'mean_workperiods'))

format = "pdf"

plt.savefig("graph_simplot." + format, 
        format=format, 
        bbox_inches="tight")
