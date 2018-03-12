#!/usr/bin/env python3

# generates the parameter files to run the simulations of the evolutionary
# reinforced threshold model



# import some libraries first
import itertools
import datetime
import shutil
import re
import sys
import pandas as pd
import numpy as np
from pathlib import Path


class RunGenerator:

    # all parameters to generate jobfiles
    # for the Carson HPC cluster
    email = "alwk201@exeter.ac.uk"
    batch_dir_prefix = "hpcbatch"
    job_file_prefix = "hpcjob"
    job_file_postfix = ".qsub"
    run_dir_prefix = "core"
    runtime_mins = 800 

    # TODO 
    # make destination directory clearer
        
    # some path necessary to compile stuff 
    ld_path = "/cm/shared/apps/gsl/gcc/1.16/lib"

    param_file_name = "params.txt"

    # initialize the class
    # note that data_dict should have columns
    # in the order of appearance in the file
    def __init__(self, 
            all_run_combinations, 
            dest_dir,
            exe):

        # store all possible combinations of parameters
        self.all_runs = all_run_combinations
   

        # store the 'destination' path where the batch 
        # folder is going to be made 
        self.dest_dir = Path(dest_dir)

        # see whether dir indeed exists
        assert(self.dest_dir.exists())

        # store the dir from which the exes have to be
        # copied
        self.exe = Path(exe)
        assert(self.exe.exists())

    # generate a batch directory with all folders
    # which contain runs and executables
    def generate_batch(self):

        # check that the dataframe has more than 0 rows
        assert(self.all_runs.shape[0] > 0)

        # folder name based on current data, so get data
        now = datetime.datetime.now()

        # set the time of creation of the parent directory
        # which will also be used as an id of the jobfiles
        idtime = now.strftime("%d_%m_%Y_%H%M%S")

        # create batch dir name
        batch_dir_name = self.batch_dir_prefix + "_" + idtime

        # make the batch folder inside the
        # destination directory
        batch_dir = self.dest_dir / batch_dir_name
        batch_dir.mkdir()

        # loop through all rows and generate the folders
        for row_number, row in self.all_runs.iterrows():

            # create new subfolder for a single run

            # make the name of the subfolder
            folder_name = self.run_dir_prefix + "_" + str(row_number)
       
            # create it
            current_folder = batch_dir / folder_name
            current_folder.mkdir()

            # write the parameter file
            self.write_parameter_file(current_folder, row_number, row)

            # copy the executable to the new directory
            shutil.copy(str(self.exe), 
                    str(current_folder / self.exe.name))

            # now make the jobfile
            self.create_jobfile(batch_dir, current_folder, row_number)

    # create the jobfile which can be submitted using qsub to run the thing.
    # batch_dir: the parent batch directory (see generate_batch())
    # run_folder: a folder within the batch dir containing a single run
    # the canonical number of that run folder
    def create_jobfile(self, 
            batch_dir, 
            run_folder, 
            run_folder_number):

        # make the name of the current job
        # this is the name that will appear in the listing generated
        # by qstat
        job_name = re.sub(str(batch_dir),"/","_slash_") + "_job_" + \
                str(run_folder_number)

        # create the start of the jobfile
        job_file_content = "#!/bin/sh\n" \
                + "#$ -N " + job_name + "\n" \
                + "#$ -S /bin/bash\n" \
                + "#$ -M " + self.email + "\n" \
                + ". /etc/profile.d/modules.sh\n" \
                + "module load shared gsl\n" \
                + "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:" \
                + self.ld_path + "\n" \
                + "cd " + str(run_folder) + "\n" \
                + "./" + self.exe.name

        # generate a filename for the qsub file
        job_file_name = self.job_file_prefix + str(run_folder_number) + \
            self.job_file_postfix

        # open the file and write the contents to it
        job_file = batch_dir / job_file_name

        job_file.write_text(data=job_file_content)


    # fill a directory with a parameter file
    # and an exe
    def write_parameter_file(self, 
            run_folder,
            run_folder_number,
            param_data):

        # first generate the contents of the file
        file_contents = ""

        for item_number, value in param_data.iteritems():
            file_contents += str(value) + "\n"

        # create the parameter file
        param_file = run_folder / self.param_file_name

        # and write to it
        param_file.write_text(data=file_contents)


    # make all possible combinations of parameter values
    # see http://pandas.pydata.org/pandas-docs/version/0.14/cookbook.html#creating-example-data 
def expand_grid(data_dict):

    # calculate total product of data frame
    rows = itertools.product(*data_dict.values())
    
    return(pd.DataFrame([row for row in rows],
        columns=data_dict.keys(),
        dtype=object))

    return(pd.DataFrame.from_records(
        rows,
        columns=data_dict.keys(),
        dtype=object))

maxtime = 1000 
# make a dictionary of all the parameters
pardict = {
        "N":[100], # number of workers / colony
        "Col": [1000], # number of colonies
        "maxtime": [maxtime], # time steps work is performed before reproduction
        "meanT1" : [ 1.0 ], # mean threshold for each task
        "meanT2" : [ 1.0 ], # mean threshold for each task
        "delta1" : [ 1.0 ], # mean threshold for each task
        "delta2" : [ 1.0 ], # mean threshold for each task
        "alpha_max_1" : [ 0.03 ], # maximum work efficiency task 1
        "alpha_max_2" : [ 0.03 ], # maximum work efficiency task 1
        "alpha_min_1" : [ 3.3e-05 ], # minimum work efficiency task 1
        "alpha_min_2" : [ 3.3e-05 ], # minimum work efficiency task 1
        "beta_1" : [ 3.3e-05 ], # stimulus decay task 1
        "beta_2" : [ 3.3e-05 ], # stimulus decay task 2 
        "p": [0.2], # quitting probability
        "mutp" : [0.01], # mutation probability
        "maxgen" : [5], # number of generations 
        "beta_fit" : [1.0], # exponent task 1 (not used)
        "gamma_fit" : [1.0], # exponent task 2 (not used)
        "recomb" : [0.5], # recombination rate
        "timecost" : [0], # duration-dependent switching cost
        "mutstep" : [0.1], # standard deviation of mutational distribution
        "initStim" : [0.3], # initial level of the stimulus
        "p_wait" : [0.1], # probability that ant has to wait c time steps before switching
        "tau" : [int(0.5 * maxtime)], # 
        "initForget" : [0.1], # 
        "initLearn" : [0.1], # 
        "step_gain_exp" : [0.2], # 
        "step_lose_exp" : [0.1], # 
        "threshold_noise" : [1.0], # 
        "K" : [0.1]
}

all_combinations = expand_grid(pardict)


all_combinations["seed"] = np.random.randint(
        low = 0, 
        high = 2147483646,
        size = all_combinations.shape[0])


# make an instance of the rungenerator class
rg = RunGenerator(
        all_run_combinations = all_combinations, 
        dest_dir="/home/bram/",
        exe="xreinforcedRT"
        )

rg.generate_batch()
