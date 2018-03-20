#!/usr/bin/env python3

# plot the division of labor simulations

import pandas as pd
import itertools
import subprocess 
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


#########################################

#           read in the data

#########################################

# get the ants' behavioral data and plot it
behav_dat = pd.read_csv("ant_beh_1.txt",
        sep=";",
        index_col=False)

#########################################

#           make figure

#########################################

# initialize the figure
fig = plt.figure(figsize=(10,18))

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


# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[0,0])

col_ids = list(behav_dat["col_id"].unique())

for col_id in col_ids:

    behav_dat_sub = behav_dat[behav_dat["col_id"] == col_id]

    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meanswitches"],
            color="grey",
            alpha=0.5)

ax.set_ylabel(r"Switches")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)


# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[1,0])

col_ids = list(behav_dat["col_id"].unique())

for col_id in col_ids:

    behav_dat_sub = behav_dat[behav_dat["col_id"] == col_id]

    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meanworkperiods"],
            color="lightblue",
            alpha=0.5)

ax.set_ylabel(r"Workperiods")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)

# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[2,0])

col_ids = list(behav_dat["col_id"].unique())

for col_id in col_ids:

    behav_dat_sub = behav_dat[behav_dat["col_id"] == col_id]

    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meanthreshold1"],
            label="_nolegend_" if col_id > 0 else r"$\theta_{1}$",
            color="blue",
            alpha=0.5)
    
    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meanthreshold2"],
            label="_nolegend_" if col_id > 0 else r"$\theta_{2}$",
            color="red",
            alpha=0.5)

ax.set_ylabel(r"$\theta_{i}$")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)

ax.legend()

# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[3,0])

for col_id in col_ids:

    behav_dat_sub = behav_dat[behav_dat["col_id"] == col_id]

    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meancountact1"],
            color="blue",
            label="_nolegend_" if col_id > 0 else r"$N_{\text{act},1}$",
            alpha=0.5)
    
    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meancountact2"],
            color="red",
            label="_nolegend_" if col_id > 0 else r"$N_{\text{act},2}$",
            alpha=0.5)

ax.set_ylabel(r"$N_{\text{act},i}$")

ax.legend()

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)

# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[4,0])

for col_id in col_ids:

    behav_dat_sub = behav_dat[behav_dat["col_id"] == col_id]

    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meanalpha1"],
            color="blue",
            label="_nolegend_" if col_id > 0 else r"$\alpha_{1}$",
            alpha=0.5)
    
    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meanalpha2"],
            color="red",
            label="_nolegend_" if col_id > 0 else r"$\alpha_{2}$",
            alpha=0.5)

ax.set_ylabel(r"$\alpha_{i}$")

ax.legend()

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)

# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[5,0])

for col_id in col_ids:

    behav_dat_sub = behav_dat[behav_dat["col_id"] == col_id]

    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meanexperiencepoints1"],
            color="blue",
            label="_nolegend_" if col_id > 0 else r"$e_{1}$",
            alpha=0.5)
    
    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["meanexperiencepoints2"],
            color="red",
            label="_nolegend_" if col_id > 0 else r"$e_{2}$",
            alpha=0.5)

ax.set_ylabel(r"Experience")

ax.legend()

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)

# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[6,0])

for col_id in col_ids:

    behav_dat_sub = behav_dat[behav_dat["col_id"] == col_id]

    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["stim1"],
            color="blue",
            label="_nolegend_" if col_id > 0 else r"$s_{1}$",
            alpha=0.5)
    
    ax.plot(
            behav_dat_sub["time"],
            behav_dat_sub["stim2"],
            color="red",
            label="_nolegend_" if col_id > 0 else r"$s_{2}$",
            alpha=0.5)

ax.set_ylabel(r"Stimulus")

ax.legend()

format = "pdf"

plt.savefig("graph_behavplot." + format, 
        format=format, 
        bbox_inches="tight")
