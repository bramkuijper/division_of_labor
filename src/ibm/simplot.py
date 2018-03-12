#!/usr/bin/env python3

# plot the division of labor simulations

import pandas as pd
import itertools
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm

plt.style.use('base')

rcParams['axes.labelsize'] = 15
rcParams['text.usetex'] = True
rcParams['font.family'] = 'sans-serif'

# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

# get colony spec data
# first read the header file
f = open("header_1.txt")
fl = f.readlines()
f.close()

header = fl[0].split("\t")

work_alloc_data = pd.read_csv("data_work_alloc_1.txt",
        sep="\t",
        header=None,
        index_col=False,
        names=header)

print(work_alloc_data.head())

sys.exit(1)

# get data on branching for learning and forgetting
branch_data = pd.read_csv(
        filepath_or_buffer=sys.argv[1], 
        sep=";", 
        header=None, # no header
        index_col=False,
        names=["minbin","maxbin","generation","trait","count"] # hence provide header names
        )

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

# generate a pivot table for forgetting
(x_forget, y_forget, forget_count) = generate_pivot(
        the_data = branch_data[branch_data["trait"]=="forget"], 
        x="generation",
        y="minbin",
        z="count"
        )


# initialize the figure
fig = plt.figure(figsize=(14,5))

# generate the grid of the graph
# see: 
widths = [ 1, 0.05 ]
heights = [ 1, 1]

# make 
gs = gridspec.GridSpec(
        nrows=len(heights),
        ncols=len(widths),
        width_ratios=widths,
        height_ratios=heights)

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




format = "pdf"

plt.savefig("graph_simplot." + format, 
        format=format, 
        bbox_inches="tight")
