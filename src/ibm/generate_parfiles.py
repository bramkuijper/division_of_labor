#!/usr/bin/env python3

# generates the parameter files to run the simulations of the evolutionary
# reinforced threshold model



# import some libraries first
import math
import itertools
import datetime
import shutil
import re
import sys
from collections import OrderedDict
import pandas as pd
import numpy as np
from pathlib import Path


class RunGenerator:

    # all parameters to generate jobfiles
    # for the Carson HPC cluster
    email = ""
    batch_dir_prefix = "hpcbatch"
    job_file_prefix = "hpcjob"
    job_file_postfix = ".qsub"
    run_dir_prefix = "core"
    runtime_mins = 800

    # after one has printed the value of a parameter
    # print a semicolon and then a label of the parameter
    print_param_key = False

    # TODO 
    # make destination directory clearer
        
    # some path necessary to compile stuff 
    ld_path = "/cm/shared/apps/gsl/2.3/lib:/cm/local/apps/gcc/7.2.0/lib"

    param_file_name = "params.txt"

    # initialize the class
    # note that data_dict should have columns
    # in the order of appearance in the file
    def __init__(self, 
            all_run_combinations, 
            dest_dir,
            exe,
            email):

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

        self.email = email

        if not email:
            print("please provide the email address")
            sys.exit(1)

    # generate a batch directory with all folders
    # which contain runs and executables
    def generate_batch(self, nrep=1):

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

        print("Creating " + str(batch_dir))


        core_number = 0

        # ok, make several replicates
        for replicate in range(0,nrep):

            # loop through all rows and generate the folders
            for rownum, row in self.all_runs.iterrows():

                # create new subfolder for a single run

                # make the name of the subfolder
                folder_name = self.run_dir_prefix + "_" + str(core_number)
           
                # create it
                current_folder = batch_dir / folder_name
                current_folder.mkdir()

                # write the parameter file
                self.write_parameter_file(current_folder, core_number, row)

                # copy the executable to the new directory
                shutil.copy(str(self.exe), 
                        str(current_folder / self.exe.name))

                # now make the jobfile
                self.create_jobfile(batch_dir, current_folder, core_number)

                core_number += 1

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

        for param_key, param_value in param_data.iteritems():
            
            if self.print_param_key:
                file_contents += str(param_value) + ";" + param_key + "\n"
            else:
                file_contents += str(param_value) + "\n"

        # create the parameter file
        # as a subdirectory of the run folder:
        param_file = run_folder / self.param_file_name

        # and write to it
        param_file.write_text(data=file_contents)


    # make all possible combinations of parameter values
    # see http://pandas.pydata.org/pandas-docs/version/0.14/cookbook.html#creating-example-data 
def expand_grid(data_dict):

    for key, value in data_dict.items():
        assert(type(value) is type([]))


    # calculate total product of data frame
    rows = itertools.product(*data_dict.values())
    
    return(pd.DataFrame([row for row in rows],
        columns=data_dict.keys(),
        dtype=object))

    return(pd.DataFrame.from_records(
        rows,
        columns=data_dict.keys(),
        dtype=object))




####################################
#  GENERATE PARAMETER FILES
####################################

# 1. this program will make a directory in your homedir
# called hpcbatch_date_time
# 2. within this directory there are folders called core_X containing
# the executable files of the Xth run of the simulation
# all executables are put in separate folders so that datafiles
# of a run are not overwritten by a simultaneous runa
# 3. the hpcbatch_date_time also contains jobfiles, called
# hpcjob_data_time_X. Each of these jobfiles contains the 
# instructions for the cluster to run one single job
# 

maxtime = 3000

# make a dictionary of all the parameters
# parameters are listed in order of appearance
# if you want to have multiple parameter values for a single
# parameter, just add them to the list, i.e., parameter1 = [ value1, value2, ..., valuen ]
pardict = OrderedDict()

maxtime = 100

pardict["N"]=[ 500] # number of workers / colony
pardict["Col"]=[ 1000] # number of colonies
pardict["maxtime"]=[maxtime] # time steps work is performed before reproduction
pardict["meanT1"]=[ 10.0 ] # mean threshold for each task
pardict["meanT2"]=[ 10.0 ] # mean threshold for each task
pardict["delta1_baseline"]=[ 1.0 ] # fixed increase in stimulus
pardict["delta2_baseline"]=[ 1.0 ] # fixed increase in stimulus
pardict["delta1"]=[ 0.0 ] # fixed increase in stimulus
pardict["delta2"]=[ 0.0 ] # fixed increase in stimulus
pardict["alfa1"]=[ 3.0 ] # maximum work efficiency task 1
pardict["alfa2"]=[ 3.0 ] # maximum work efficiency task 1
pardict["beta1"]=[ 0.0 ] # maximum work efficiency task 1
pardict["beta2"]=[ 0.0 ] # maximum work efficiency task 1
pardict["exp_task_1"]=[0.5] # exponent task 1 
pardict["exp_task_2"]=[0.5] # exponent task 2 
pardict["p"]=[ 0.2] # quitting probability
pardict["mutp"]=[0.1] # mutation probability
pardict["mutstd"]=[0.1] # mutation probability
pardict["recomb"]=[0,0.5] # mutation probability
pardict["maxgen"]=[10000] # number of generations 

# number of timesteps a worker has to wait 
# before engaging in another task (this is a cost of switching)
pardict["timecost"]=[10] 

# initial value of the stimulus
pardict["initStim"]=[0] 
pardict["p_wait"]=[1.0] 
pardict["tau"]=[int(maxtime/2)] 
pardict["A"]=[0.0] #Deterministic factor
pardict["B"]=[0.0] #Stochastic factor  
pardict["genspercycle"]=[10, 20, 30, 40, 50] #Generations per environmental cycle      
#pardict["randommax"]=[0, 0.001, 0.01, 0.1, 1, 10] #Maximum value of positive random number
pardict["randommax"]=[0] #Maximum value of positive random number


# make all parameter combinations
# this can be left alone
all_combinations = expand_grid(pardict)


# add a column with random numbers representing the seed
# this can be left alone
all_combinations["seed"] = np.random.randint(
        low = 0, 
        high = 2147483646,
        size = all_combinations.shape[0])

# make an instance of the rungenerator class
# change stuff here
rg = RunGenerator(
        all_run_combinations = all_combinations, 
        dest_dir=str(Path.home()), # put hpcbatch in home directory
        exe="xfixed_response", # SET THE EXECUTABLE HERE
        email="ngt206@exeter.ac.uk"
        )

# generate the batch with the number of replicates per batch
# if you just want a single replicate for each parameter combination
# nrep=1 should do it
rg.generate_batch(nrep=1)
