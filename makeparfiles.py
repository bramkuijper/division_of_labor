#!/usr/bin/env python
#python script to make parameter and job files

import sys
import os
import random
import shutil
import datetime

myprog = "stepsize"

workers = 100  
col = 1000
steps = 3000
thresh1 = 10.0
thresh2 = 10.0 
delta1 = 1.0
delta2 = 1.0
alfamax1 = 10.0
alfamax2 = 10.0
alfamin1 = 3.0
alfamin2 = 3.0
beta1 = 0.0
beta2 = 0.0
quitp = 1 
mutp = 0.01
gen = 10000 
beta_fit = 0.5
#gama_fit = 0.5
#seed = 
recomb = 1 
timecost = [0, 1, 2, 4, 6,8,10]
mutstep = 0.01
initstim = 10 
p_wait = 1
tau = 10 
initForget = 0
initLearn = 0
step_gain_exp = [0, 1, 5, 10, 15] 
step_lose_exp =  [0, 1, 5, 10, 15]
K = 0.05 

rep = range(0,10,1)

random.seed()


# executable path
executable_path  = "/home/uccoaku/work_ana/anas_stuff"

# prefix of the directory in which jobs are run
parent_dir_prefix = "hpcbatch_ana"

# prefix of the jobfile that is submitted via qsub
jobfile_prefix = "hpcjob"

# name identifier of the job
job_name = "division_labor"

# project identifier
project_id_str = "mitochondria"

# memory
memory_str = "10M"

# email
email = "a.kuijper@ucl.ac.uk"

# working directory
work_dir = "/home/uccoaku/Scratch"

walltime_hours = "20"
walltime_minutes = "00"


def make_parameter_files():
    i=0
    dir =[]

    # get the current date and time
    now = datetime.datetime.now()

    main_dir = work_dir + "/" + parent_dir_prefix + "_" + now.strftime("%d_%m_%Y_%H%M%S")

    os.mkdir(main_dir)
    os.chdir(parentdir)

    for b in beta_fit:
        for c in timecost:
            for gain in step_gain_exp:
                for lose in step_lose_exp:
                    for r in rep:

                        seed = random.randint(1,5000)
                        my_pars=[]
                        my_job=[]
                        gama_fit = 1.0 - b
                        dirName = main_dir + "/" + "reinfRT_gain_" + str(gain) +"_lose_"+str(lose)+"_cost_"+str(c) + "_beta_" + str(b)+ "_rep_" +str(r)
                        dir.append(dirName)
                       
                        if not os.path.isdir("./" + dirName + "/"):
                            os.mkdir("./" + dirName + "/")

                        #create parameter list
                        my_pars.append(str(workers)+"\n") 
                        my_pars.append(str(col)+"\n")
                        my_pars.append(str(steps)+"\n")
                        my_pars.append(str(thresh1)+"\n")
                        my_pars.append(str(thresh2)+"\n") 
                        my_pars.append(str(delta1)+"\n") 
                        my_pars.append( str(delta2)+"\n")
                        my_pars.append( str(alfamax1)+"\n")
                        my_pars.append( str(alfamax2)+"\n")
                        my_pars.append( str(alfamin1)+"\n")
                        my_pars.append( str(alfamin2)+"\n")
                        my_pars.append( str(beta1)+"\n")
                        my_pars.append( str(beta2)+"\n")
                        my_pars.append( str(quitp)+"\n")
                        my_pars.append( str(mutp)+"\n")
                        my_pars.append( str(gen)+"\n")
                        my_pars.append( str(b)+"\n")
                        my_pars.append( str(gama_fit)+"\n")
                        my_pars.append( str(seed)+"\n")
                        my_pars.append( str(recomb)+"\n")
                        my_pars.append( str(c)+"\n")
                        my_pars.append( str(mutstep)+"\n")
                        my_pars.append( str(initstim)+"\n")
                        my_pars.append( str(p_wait)+"\n")
                        my_pars.append( str(tau)+"\n")
                        my_pars.append( str(initForget)+"\n")
                        my_pars.append( str(initLearn)+"\n")
                        my_pars.append( str(gain)+"\n")
                        my_pars.append( str(lose)+"\n")
                        my_pars.append( str(K)+"\n")
                        
                        filename = dirName + "/params.txt"
                        fileobj = open(filename, 'w')
                        fileobj.writelines(my_pars) 
                        fileobj.close()
                       
                        shutil.copy(myprog, dirName)            

                        # generate the jobfile's name
                        jobname= dirname + "/" + jobfile_prefix + "_" + "reinfRT_gain_" + str(gain) +"_lose_"+str(lose)+"_cost_"+str(c) + "_beta_" + str(b)+ "_rep_" +str(r) + "_" + ".qsub"

                        # we also need this jobfile name if we want to use jobtracker
                        # but then in a different format, with slashes replaced by pluses
                        jobfilename_sub = re.sub("\/","+",jobname)
                        
                        # now generate the contents of the jobfile
                        my_job_str =  "#!/bin/bash -l\n \
                                        #$ -S /bin/bash\n \
                                        #$ -l h_rt=" + str(walltime_hours) + ":" + str(walltime_minutes) + ":00\n \
                                        #$ -l mem=" + memory_str + "\n \
                                        #$ -l thr=1\n \
                                        #$ -N " + jobfilename_sub + "\n \
                                        #$ -P " + project_id_str + "\n \
                                        #$ -wd " + dirName + "\n \
                                        cd $TMPDIR \n" \
                                        + executable_path + "\n \
                                        find . -maxdepth 1 -type f -iname \"*.txt\" -exec mv {} " \
                                        + dirName + "/. \;" + "\n"
                        
                        
                        jobfile = open(jobname, "w")                
                        jobfile.write(my_job_str)
                        jobfile.close()
                         
                        i+=1
make_parameter_files()


