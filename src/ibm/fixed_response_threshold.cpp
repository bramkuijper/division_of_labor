
// with stimulus update per ant, or per timestep (see define)

//---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/stat.h>
#include <sstream>

//#define DEBUG
//#define SIMULTANEOUS_UPDATE
//#define STOPCODE
//---------------------------------------------------------------------------

using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *rng_global; // gnu scientific rng 


struct Params
{
    int N; //number of workers
    int Col; // number of colonies
    int maxtime; // time steps
    double p;  // the probability that a working ant quits working
    int tasks; // number of tasks
    int maxgen; // number of generations 
    // the seed for the random number generator
    int seed;

    // the recombination rate between the threshold loci
    // coding for the different tasks
    // (scaled between 0 [no recombination]
    // to 0.5 [full recombination]
    double recomb;
    
    // the number of timesteps ants cannot work when switching
    // between different tasks
    int timecost;

    double mutp;//mutation probability
    double mutstd; // standard deviation of the mutational distribution 

    double initStim; // initial stimulus value
    double p_wait; // probability that ant has to wait c time steps before switching
    int tau; // time step from which fitness is counted

    vector<double> meanT; // mean thresholds
    vector<double> delta; // variable to keep current stimulus increase level
    vector<double> deltabaseline; // baseline increase in delta

    // efficiency with which workers perform 
    // tasks (see Bonabeau et al 1996 eq. 3)
    vector<double> alfa; 

    // the decay parameter in stimulus 
    // increase, see UpdateStim()
    vector<double> beta; 

    // exponents that weigh how strongly each task impinges on fitness
    // these are the beta and 1-beta parameters 
    // in eq (3) of Duarte et al 2012 Behav Ecol Sociobiol
    vector<double> fitness_weights; 


    // stochsine
    double A; //Deterministic factor
	double B; //Stochastic factor
	int genspercycle; //Generations per environmental cycle
	int randommax; //Maximum value of positive random number
	int gensdone; //Generations completed
	int stepsdone; //Timesteps completed in current generation
  

    // function to initialize the parameters from the parameter
    // file
    istream & InitParams(istream & inp);
}; // en strut params

struct Ant
{
    // genome
    vector <double> threshold;  // thresholds of this ant
    vector <int> countacts;   // counter of acts done by this ant
    vector < bool > want_task; // TODO
    int last_act; // the ant's last act (needed to assess whether it switches or not)
    int curr_act; // the ant's current act
    int switches; // number of transitions to a different task
    int workperiods; // number of working periods 
    double F; // specialization value
    bool mated; // whether queen has been mated already
    int count_time; // counter of timesteps to switch task
};

typedef vector <Ant>  Workers;
typedef vector <Ant> Sexuals;

struct Colony
{
    Workers MyAnts; // stack of workers
    Ant male, queen; // queen and her male
    int ID; // unique ID of the colony

    vector<double> stim; // stimulus level at time t for each task
    vector<double> newstim; // stimulus level at time t+1 for each task
    vector<double> workfor;  // number acts * eff each time step
    vector <int> numacts; // number of acts performed per task

    double idle; // proportion workers that _never_ worked in the simulation 
    double inactive; // proportion workers that were idle each time step
    vector <double>last_half_acts; //number of acts performed in the last half of simulation

    // product of number of acts for each task, 
    // every time step (of half of simulation)
    double work_fitness; 

    // TODO
    double fitness;

    double diff_fit;// fitness difference to minimal fitness
    double rel_fit; // fitness relative to whole population
    double cum_fit; //cumulative fitness

    // statistics
    vector<double> HighF, LowF, CategF;
    vector<double> HighT, LowT, CategT1, CategT2;
    vector<double> mean_work_alloc;
    double mean_F;
    double var_F;
    double mean_F_franjo;
    double var_F_franjo;
    int num_offspring;
    double mean_switches;
    double var_switches;
    double mean_workperiods;
    double var_workperiods;
};

// specify the population
typedef vector < Colony > Population;

// vector containing selected colonies that will 
// reproduce
vector <int> parentCol;
double sum_Fit;
Sexuals mySexuals;


int simstart_generation;
int simpart;


// A for loop/stream that reads in the parameter file, repeating for the number of tasks.
istream & Params::InitParams(istream & in)
{
    // set the number of tasks to 2
    tasks = 2;

    // allocate space in vectors
    meanT.reserve(tasks);
    delta.reserve(tasks);
    alfa.reserve(tasks);
    beta.reserve(tasks);
    fitness_weights.reserve(tasks);

    // read in the first three parameters
    in >> N >>
        Col >>
        maxtime; 

    // temporary variable to process current input
    double tmp;

    // read in the initial mean thresholds for each stimulus
    for (int i = 0; i < tasks; ++i) 
    {
        in >> tmp; 
        meanT.push_back(tmp);
    }

    // read in the deltas
    // this will change with stochsine TODO 
    for (int i = 0; i < tasks; ++i) 
    {    
        in >> tmp ; 
        deltabaseline.push_back(tmp);
    }
    
    // read in the deltas
    // this will change with stochsine TODO 
    for (int i = 0; i < tasks; ++i) 
    {    
        in >> tmp ; 
        delta.push_back(tmp);
    }

    // read in the work efficiency parameter (see eq. 3 in Bonabeau et al 1996)
    for (int i = 0; i < tasks; ++i) 
    {
        in >> tmp; 
        alfa.push_back(tmp);
    }

    // read in the basic growth rate of the stimulus, see UpdateStim
    // Basically, we allow for s(t+1) = (1-beta)*s(t) + delta - alpha * f(Nworkers);
    // whereas Bonabeau et al assume s(t+1) = s(t) + delta - alpha * f(Nworkers)
    for (int i = 0; i < tasks; ++i) 
    {
        in >> tmp;
        beta.push_back(tmp);
    }

    // read in the fitness exponents used for each task
    // see Duarte et al 2012 Behav Ecol Sociobiol eq (3)
    for (int i = 0; i < tasks; ++i) 
    {
        in >> tmp;
        fitness_weights.push_back(tmp);
    }

    in >> p >>
        mutp >>
        mutstd >>
        recomb >>
        maxgen >>
        timecost>>
	initStim >>
	p_wait >>
	tau >>

    //stochsine			 
    A >> //Deterministic factor
    B >> //Stochastic factor
    genspercycle >> //Generations per environmental cycle
    randommax >> //Maximum value of positive random number
    seed;

    return(in);
}

// function to check if a file exists already
bool FileExists(string strFilename) 
{
    struct stat stFileInfo;
    int intStat;

    // Attempt to get the file attributes
    intStat = stat(strFilename.c_str(),&stFileInfo);

    if (intStat == 0) 
    {
        // We were able to get the file attributes
        // so the file obviously exists.
        return(true);
    } 

    return(false);
}

//========================================================================================
// Function that initializes founders from (previous) data 
void StartFromLast(istream &in, Params &Par, Population &Pop)
{
    // temporary variable to process current input
	double tmp;
	
	in >> simpart 
	   >> simstart_generation; 

	for (unsigned int i = 0; i < Pop.size(); ++i)
    {
        for (int task = 0; task < Par.tasks; ++task) 
        {
            in >> tmp;
            Pop[i].male.threshold[task] = tmp;
        }

        for (int task = 0; task < Par.tasks; ++task) 
        {
            in >> tmp;
            Pop[i].queen.threshold[task] = tmp;
        }
    }
}

// Diagnostic function: prints the results of the par object 
// to verify that reading in the parameters was performed correctly.
void ShowParams(Params & Par)
{
    cout <<"Workers " << Par.N << endl;
    cout << "Colonies " << Par.Col << endl;
    cout << "Timesteps " << Par.maxtime << endl;
 
     for (int task = 0; task < Par.tasks; ++task)
     {
        cout << "Initial T" << task << "\t" << Par.meanT[task] << endl;
     }
 
     for (int task = 0; task < Par.tasks; ++task)
     {
        cout << "Delta " << task << "\t"  << Par.delta[task] << endl;
     }
     
     for (int task = 0; task < Par.tasks; ++task)
     {
        cout << "Delta" << task << "baseline\t"  << Par.deltabaseline[task] << endl;
     }
 
     for (int task = 0; task < Par.tasks; ++task)
     {
        cout << "Efficiency " << task << "\t" << Par.alfa[task] << endl;
     }
 
    for (int task = 0; task < Par.tasks; ++task)
    {
        cout << "Decay " << task << "\t" << Par.beta[task] << endl;
    }
    
    for (int task = 0; task < Par.tasks; ++task)
    {
        cout << "Fitness weight " << task << "\t" << Par.fitness_weights[task] << endl;
    }
 
    cout << "prob quit" << Par.p << endl;
    cout << "Task number " << Par.tasks << endl;
    cout << "mut prob " << Par.mutp << endl;
    cout << "mut std " << Par.mutstd << endl;
    cout << "recombination " << Par.recomb << endl;
    cout << "Max gen " <<  Par.maxgen << endl;
    cout << "seed " << Par.seed << endl;
    cout << "timecost " << Par.timecost << endl;
    cout << "tau " << Par.tau << endl;
    cout << "Deterministic Factor A " << Par.A << endl;
    cout << "Stochastic Factor B " << Par.B << endl;
    cout << "Generations per Cycle " << Par.genspercycle << endl;
    cout << "Maximum Random Number " << Par.randommax << endl;
}
//-----------------------------------------------------------------------------



// Initialises the founding generation by mating a male and a queen, with attributes randomly selected from a normal distribution
void InitFounders(Population &Pop, Params &Par)
{
    Pop.resize(Par.Col);

    // initialize thresholds of males and queens
	for (unsigned int i = 0; i < Pop.size(); ++i)
    {
        Pop[i].male.threshold.resize(Par.tasks);
        Pop[i].queen.threshold.resize(Par.tasks);

        for (int task = 0; task < Par.tasks; ++task)
        {
            // TODO
            // why is this pop.size/2??
            if (i < Pop.size()/2)
            {
                Pop[i].male.threshold[task] = Par.meanT[task];
                Pop[i].queen.threshold[task] = Par.meanT[task];
		    }

            Pop[i].male.mated =true;
            Pop[i].queen.mated = true;
		}
	}
}// end InitFounders

// mutation of a single threshold allele
double Mutate(double val, double mu, double mustd)
{
    if (gsl_rng_uniform(rng_global) < mu)
    {
        val += gsl_ran_gaussian(rng_global, mustd);
    }

    // thresholds cannot be <0 
    if (val < 0)
    {
        val = 0;
    }

    return(val);
}


// Defines inheritance, producing worker 
// threshold genotypes from their parents, 
// including mutation
void Inherit(Ant &Daughter, Ant &Mom, Ant &Dad, Params &Par)
{
    // start with inheritance from mom (true) or from dad (false)
    bool inherit_from_mom = gsl_rng_uniform(rng_global) < 0.5;

    double allelic_value;

    // now start to inherit all the thresholds
    for (int task = 0; task < Par.tasks; ++task)
    {
        // when beyond the first gene locus, 
        // check whether we have recombined
        if (task > 0)
        {
            // ok recombination happened, hence change parent 
            // from which the next allele will originate
            if (gsl_rng_uniform(rng_global) < Par.recomb)
            {
                inherit_from_mom = !inherit_from_mom;
            }
        }

        // inherit the allele from the designated parent
        allelic_value = inherit_from_mom ? Mom.threshold[task] : Dad.threshold[task];

        // mutate it
        allelic_value = Mutate(allelic_value, Par.mutp, Par.mutstd);

        // assign it to daughter
        Daughter.threshold[task] = allelic_value;
    }
} // end of Inherit


// Initialises the ant workers, as generated by the Inherit function
// Prints their threshold variables, and puts them to work on a task depending on their inherited thresholds
void InitAnts(Ant & myAnt, Params & Par, Colony & myCol)
{
    // allocate space for the thresholds
    myAnt.threshold.resize(Par.tasks);

    if (Par.maxgen > 1)
    {
        // inherit from mom and dad
        Inherit(myAnt, myCol.queen, myCol.male, Par);
    }
    else 
    {
        // or just assing thresholds from parameters
        for (int task = 0; task < Par.tasks; ++task)
        {
            myAnt.threshold[task] = Par.meanT[task];
        }
    }

    myAnt.countacts.resize(Par.tasks);
    myAnt.want_task.resize(Par.tasks);

    // initialize counters and whether ant wants to perform a task
    for (int task = 0; task < Par.tasks; ++task)
    {
        myAnt.countacts[task] = 0;
        myAnt.want_task[task] = false;
    }
    //other stuff
    myAnt.last_act = Par.tasks; // 0 = task 1, 1 = task 2, etc, n = idle
    myAnt.curr_act = Par.tasks; // 0 = task 1, 1 = task 2, etc, n = idle
    myAnt.switches = 0; // swich counter of switches to 0
    myAnt.workperiods=0;
    myAnt.F = 0;
    myAnt.mated = false;
    myAnt.count_time=0;
}


// Initilializes all colonies at the start of each evolutionary generation
void Init(Population & Pop, Params & Par)
{
    for (unsigned int colony_i = 0; colony_i < Pop.size(); ++colony_i)
    {
        Pop[colony_i].MyAnts.resize(Par.N);
        Pop[colony_i].ID = colony_i;
        Pop[colony_i].fitness = 0;
        Pop[colony_i].rel_fit = 0;
        Pop[colony_i].cum_fit = 0;
        Pop[colony_i].diff_fit = 0;
        Pop[colony_i].idle = 0;
        Pop[colony_i].inactive = 0;
        Pop[colony_i].mean_F = 10; // specialization value
        Pop[colony_i].var_F=0;
        Pop[colony_i].mean_F_franjo = 10;
        Pop[colony_i].var_F_franjo = 0;
        Pop[colony_i].mean_switches = 0;
        Pop[colony_i].var_switches = 0;
        Pop[colony_i].mean_workperiods=0;
        Pop[colony_i].var_workperiods=0;
        Pop[colony_i].work_fitness = 0;

        // erase any existing workforce numbers
        Pop[colony_i].workfor.clear();
        Pop[colony_i].last_half_acts.clear();
        Pop[colony_i].numacts.clear();
        Pop[colony_i].stim.clear();
        Pop[colony_i].newstim.clear();
        Pop[colony_i].mean_work_alloc.clear();

        // allocate space for them
        Pop[colony_i].workfor.reserve(Par.tasks);
        Pop[colony_i].last_half_acts.reserve(Par.tasks);
        Pop[colony_i].numacts.reserve(Par.tasks);
        Pop[colony_i].stim.reserve(Par.tasks);
        Pop[colony_i].newstim.reserve(Par.tasks);
        Pop[colony_i].mean_work_alloc.reserve(Par.tasks);

        for (int job = 0; job < Par.tasks; ++job)
        {
            Pop[colony_i].workfor.push_back(0);
            Pop[colony_i].last_half_acts.push_back(0);
            Pop[colony_i].numacts.push_back(0);
            Pop[colony_i].stim.push_back(Par.initStim);
            Pop[colony_i].newstim.push_back(0);
            Pop[colony_i].mean_work_alloc.push_back(0);
        }
    
        for (unsigned int j = 0; j < Pop[colony_i].MyAnts.size(); ++j)
        {
            InitAnts(Pop[colony_i].MyAnts[j], Par, Pop[colony_i]);
        }
    } // end of for colony_i
} // end of Init()

// diagnostic function to show all individual ants of a colony
void ShowAnts(Colony & anyCol)
{
    cout << "Current values of:" << endl; 

    for (unsigned int ant = 0; ant < anyCol.MyAnts.size(); ++ant)
	{
	    cout << "ant " << ant << endl;
	    cout << "F " << anyCol.MyAnts[ant].F << endl; // specialization value
	    cout << "switches " << anyCol.MyAnts[ant].switches << endl;
	    cout << "workperiods " << anyCol.MyAnts[ant].workperiods << endl;
	    cout << "count acts " << anyCol.MyAnts[ant].countacts[0] 
            << "\t" << anyCol.MyAnts[ant].countacts[1]  << endl;
	}
}

void ShowColony(Colony & anyCol)
{
    cout << "Current values of:" << endl; 
    cout << "fitness " << anyCol.fitness << endl;
    cout << "relative fitness " << anyCol.rel_fit << endl;
    cout << "cumulative fitness " << anyCol.cum_fit << endl;
    cout << "diff fit " << anyCol.diff_fit << endl;
    cout << "idle " << anyCol.idle << endl;
    cout << "inactive " << anyCol.inactive << endl;
    cout << "mean_F " << anyCol.mean_F << endl; // specialization value
    cout << "var_F" << anyCol.var_F << endl;
    cout << "mean_franjo" << anyCol.mean_F_franjo << endl;
    cout << "var_franjo " << anyCol.var_F_franjo << endl;
    cout << "mean_switches " << anyCol.mean_switches << endl;
    cout << "var_switches " <<anyCol.var_switches << endl;
    cout << "workperiods " << anyCol.mean_workperiods << endl;
    cout << "var workperiods " <<  anyCol.var_workperiods << endl;
	cout << "workfor 1 " <<  anyCol.workfor[0]/3 << endl;
	cout << "workfor 2" <<  anyCol.workfor[1]/3 << endl;
	ShowAnts(anyCol);
}

// Updates the stimulus per ant for each timestep.
// Stimulus for a task decreased based of how much work put towards it.
//
// funtion is only called when there is no simultaneous updating of 
// the stimulus levels
void UpdateStimPerAnt(Params & Par, Colony & anyCol, Ant & anyAnt, int task)
{
    anyCol.workfor[task] += Par.alfa[task];   //update the amount of work done (fitness) 
    anyCol.stim[task] -= (Par.alfa[task]/Par.N); //update stimulus     
    
    if (anyCol.stim[task] < 0)
    {
        anyCol.stim[task] = 0;
    }
}


// Checks whether an ant will cease task performance and becomes idle
// Based on random draws of an ant's 
// innate quitting probability, unrelated to any other varaiable
void QuitTask(Colony & anyCol, Ant & anyAnt, int job, Params & Par)
{
    // ant quits
    if (gsl_rng_uniform(rng_global) < Par.p)
    {
        // ant does not want to do current task
        anyAnt.want_task[anyAnt.curr_act] = false;
        anyAnt.count_time = 0; // reset time to zero, she may choose the same or another task next
        anyAnt.curr_act = Par.tasks; // set to no task currently
    }
    else // ant works on
    {
        // if no simultaneous update, update the stimulus levels for
        // this ant
#ifndef SIMULTANEOUS_UPDATE
        UpdateStimPerAnt(Par, anyCol, anyAnt, job);
#endif
    }

}

// Defines response probability to certain stimulus, influenced by threshold genotype
double RPfunction (double t, double s) 
{
    double RP;

    if (t > 0.00) 
    {
        // positive threshold, calculate 
        // response according to Bonabeau et al eq 1
        RP = s*s / ((s*s)+(t*t));
    }
    else if (s > 0.00)
    {
        RP = 1.0; // if threshold is zero, then probability of acting is 1.  
    }
    else
    {
        RP = 0.0;
    }
    
    return(RP);
}
//-------------------------------------------------------------------------------

// take account of which task an ant is currently engaged
// in (and in absence of simultaneous updating)
void DoTask(Params Par, Colony & anyCol, Ant & anyAnt,int job)
{
    // updating her current act for the task she's doing
    anyAnt.curr_act = job; 

    // increasing her workperiods 
    ++anyAnt.workperiods;

#ifndef SIMULTANEOUS_UPDATE 
    //updating stimulus immediately 
    UpdateStimPerAnt(Par, anyCol, anyAnt, job);
#endif
}
//end DoTask
//

// check which stimuli (when noise added) exceed the thresholds
void WantTask(Params Par, Colony & anyCol, Ant & anyAnt)
{
    // variable that stores which tasks
    // have stimulus levels that exceed the threshold
    vector<int>counter;
    counter.reserve(Par.tasks); 

    for (int task_i = 0; task_i < Par.tasks; ++task_i)
	{
        // add random noise to both threshold and stimulus
        double stim_noise = anyCol.stim[task_i] + 
            gsl_ran_gaussian(rng_global,1.0);

        double t_noise =  anyAnt.threshold[task_i] + 
            gsl_ran_gaussian(rng_global,1.0);

        if (stim_noise < 0)
        {
            stim_noise = 0;
        }

        if (t_noise < 0)
        {
            t_noise = 0;
        }

        // check which tasks have stimulus levels beyond
        if (stim_noise >= t_noise && stim_noise > 0) 
        {
            // the threshold
            counter.push_back(task_i);  
        } 
        else 
        {
            // if her threshold not high enough, then quit with wanting task
            anyAnt.want_task[task_i] = false; 
        }
	}

    // if more than one task is above threshold, 
    // do random task among those
    if (counter.size() > 1) 
	{
        // select a random job
		int job = gsl_rng_uniform_int(rng_global, counter.size());
		anyAnt.want_task[counter[job]] = true;
	}
    else if (!counter.empty())
    {
        // otherwise only want the firs (and only) task
        anyAnt.want_task[counter[0]] = true;
    }

    // if stimulus levels exceed all the thresholds
    // ant will not want do any job
}
//end WantTask


// evaluate whether ant can switch to task it wants to do
void EvalTaskSwitch (Params & Par, Colony & anyCol, Ant & anyAnt, int myjob)
{
    // if current job is the last job
    // that this ant did, there is no switching
    // let ant just get on with the job
    if (myjob == anyAnt.last_act) 
    {
        // set job stats and update stimuli
        DoTask(Par, anyCol, anyAnt, myjob);
    }
    else // ant performed a different task relative to what it wants to do now
    {
        // find out whether ant cannot switch 
        // to a different ask but has to wait
        if (gsl_rng_uniform(rng_global) <= Par.p_wait
                && anyAnt.count_time < Par.timecost)    
        {
            anyAnt.curr_act = Par.tasks; // stays idle for as long as count_time<timecost    
            ++anyAnt.count_time; // increase the time counter as ant is doing nothing
        }
        else // ok ant can switch tasks now
        {
            // set job stats and update stimuli
            DoTask(Par, anyCol, anyAnt, myjob);
        }
    }
}
//end EvalTaskSwitch


// Causes an ant to perform one of the tasks, using responce probability
// Involves the stimulus per ant and threshold, and defines when an ant might switch tasks
void TaskChoice(Params & Par, Colony & anyCol, Ant & anyAnt)
{ 
    // check if ant previously did not want to do neither task 
    if (!anyAnt.want_task[0] && !anyAnt.want_task[1])
    {
        // assess whether the ant still does not want to do
        // any tasks
        WantTask(Par, anyCol, anyAnt);
    }

    // go through the tasks and see whether they want to be done
    for (int task = 0; task < Par.tasks; ++task)
    {
        // ok ant wants to do a task
        // (it can only want to do one task)
        if (anyAnt.want_task[task])
        {
            // ant currently doing nothing
            if (anyAnt.last_act >= Par.tasks) 
            {
                // perform the task
                DoTask(Par, anyCol, anyAnt, task);
            }
            else // ant currently engaged in a task
            {
                EvalTaskSwitch(Par, anyCol, anyAnt, task); 
            }
        }    
        else 
        {
            // if she doesn't want the task, just remain idle
            anyAnt.curr_act = Par.tasks; 
        }
    } // end for task_i
} // end of TaskChoice()


// update ants and number of switches
void UpdateAnts(Population & Pop, Params & Par)
{
    int current_act;

    // loop through colonies
    for (unsigned int colony_i = 0; colony_i < Pop.size(); ++colony_i)
    {
        Pop[colony_i].inactive = 0;

        for (int task = 0; task < Par.tasks; ++task)
        {
            // reset statistics on work for tasks to 0
            Pop[colony_i].workfor[task] = 0; 
        }
        // go through individual ants of the colony
        //
        // let active ants potentially quit
        // let idle ants potentially find work

        // first shuffle vectors
        random_shuffle(Pop[colony_i].MyAnts.begin(), Pop[colony_i].MyAnts.end());
     
        for (unsigned int ant_i = 0; 
                ant_i < Pop[colony_i].MyAnts.size(); ++ant_i)  
        {
            // check whether ant is currently active
            if (Pop[colony_i].MyAnts[ant_i].curr_act < Par.tasks)
            {
                Pop[colony_i].MyAnts[ant_i].last_act = 
                    Pop[colony_i].MyAnts[ant_i].curr_act; //record last act

                // check whether ant quits
                QuitTask(Pop[colony_i], 
                        Pop[colony_i].MyAnts[ant_i], 
                        Pop[colony_i].MyAnts[ant_i].curr_act, 
                        Par); 
            }

            // ant is currently inactive
            if (Pop[colony_i].MyAnts[ant_i].curr_act >= Par.tasks)
            {
                //if inactive, choose a task 
                TaskChoice(Par, Pop[colony_i], Pop[colony_i].MyAnts[ant_i]); 
            }

            //if ant (still nor just now) active 
            // update counters
            if (Pop[colony_i].MyAnts[ant_i].curr_act < Par.tasks)
            {
                current_act = Pop[colony_i].MyAnts[ant_i].curr_act;
                // update the number of acts done
                ++Pop[colony_i].numacts[current_act];
                ++Pop[colony_i].MyAnts[ant_i].countacts[current_act];
           
                // update number of switches
                if (Pop[colony_i].MyAnts[ant_i].last_act < Par.tasks 
                        && Pop[colony_i].MyAnts[ant_i].last_act != 
                        Pop[colony_i].MyAnts[ant_i].curr_act)
                {
                    ++Pop[colony_i].MyAnts[ant_i].switches;
                }
            }
            else
            {
                ++Pop[colony_i].inactive;
            }
        }// ant for ant_i 

        // update proportion of inactive workers 
        Pop[colony_i].inactive /= Pop[colony_i].MyAnts.size(); 

    } // end for colony_i
}  // end of UpdateAnts()
//------------------------------------------------------------------------------

//stochsine
//Creates value for delta with a stochastic sine wave (Botero et al. 2015)
void Stochsine(Params & Par)
{
    for (int task_i = 0; task_i < Par.tasks; ++task_i)
    {
        Par.delta[task_i] =  //plus one for baseline stimulus increase
            Par.deltabaseline[task_i] + (Par.A  * sin((2 *

            //pi
            M_PI *

            //Calculate cumulative timesteps
            ((Par.maxtime * Par.gensdone) + Par.stepsdone)

            ) / Par.maxtime* Par.genspercycle))
            + (Par.B *

            //Random number between 0 and randdommax
            gsl_rng_uniform_pos(rng_global) * Par.randommax);
    }
}



// Increases the stimulus by delta
// Decreases the stimulus depending on the amount of work done towards a task
void UpdateStim(Population & Pop, Params & Par)   
{
	//stochsine
	Stochsine(Par);

    // go through all colonies
    for (unsigned int colony_i = 0; colony_i < Pop.size(); ++colony_i)
    {

        for (int task=0; task < Par.tasks; ++task)
        {
            // update the value of the stimulus as in 
            // in eq. (3) of Bonabeau, with the difference that 
            // s(t) is multiplied by decay parameter 1-beta
            Pop[colony_i].newstim[task] = (1.0 - Par.beta[task]) * Pop[colony_i].stim[task] 
                + Par.delta[task];
            
            
#ifdef SIMULTANEOUS_UPDATE
            // in case of simultaneous update subtract all the work done
            // from the stimulus dynamic in one go
            // otherwise this will be done per ant later on
            Pop[colony_i].newstim[task] -= (Pop[colony_i].workfor[task]/Par.N); 
#endif

            // stimulus cannot be negative
            if (Pop[colony_i].stim[task] < 0)
            {
                Pop[colony_i].stim[task] = 0;
            }

            // update the stimulus
            Pop[colony_i].stim[task] = Pop[colony_i].newstim[task];
        } // end for task
    } // end for colony_i
} // end UpdateStim()
//------------------------------------------------------------------------------

void Calc_F(Population & Pop, Params & Par) // calculate specialization 
{
    double C;

    for (unsigned int colony_i = 0;  colony_i < Pop.size(); ++colony_i)
    {
        Pop[colony_i].mean_F = 0; // F varies between -1 and 1
        Pop[colony_i].mean_F_franjo = 0; // F_franjo varies between 0 and 1 
        Pop[colony_i].mean_switches = 0; 
        Pop[colony_i].mean_workperiods = 0; 

        double F_franjo;
        double mean_F_franjo = 0;

        double sumsquares_F = 0;
        double sumsquares_F_franjo = 0;

        double sumsquares_switches = 0;
        double sumsquares_workperiods = 0;

        // keep track of the total number of ants
        // who have worked at least once
        size_t n_ants_active = 0;

        for (unsigned int ant_i = 0; ant_i < Pop[colony_i].MyAnts.size(); ++ant_i)
        {
            assert(Pop[colony_i].MyAnts[ant_i].workperiods <= Par.maxtime);
            C = 0;

            // calculate average number of workperiods
            Pop[colony_i].mean_workperiods += Pop[colony_i].MyAnts[ant_i].workperiods;
            
            // calculate sum of squares for variance
            sumsquares_workperiods += Pop[colony_i].MyAnts[ant_i].workperiods * 
                Pop[colony_i].MyAnts[ant_i].workperiods;

            // if this ant has ever worked take statistics on switches
            if (Pop[colony_i].MyAnts[ant_i].workperiods > 0)
            {
                // calculate mean switches and 
                // workperiods for book-keeping
                Pop[colony_i].mean_switches += Pop[colony_i].MyAnts[ant_i].switches;

                // calculate sum of squares for variance
                sumsquares_switches += Pop[colony_i].MyAnts[ant_i].switches * 
                    Pop[colony_i].MyAnts[ant_i].switches;
      

                //  C is frequency of switching between tasks 
                //  which is the number of switches divided by the total number of
                //  possible switches, which is the number of workperiods minus 1
                //  (as one cannot switch anymore during the final workperiod)
                C = double(Pop[colony_i].MyAnts[ant_i].switches) / 
                    (Pop[colony_i].MyAnts[ant_i].workperiods - 1.0);
     
                // F is between -1 and 1
                Pop[colony_i].MyAnts[ant_i].F = 1 - 2*C;
                // F_franjo is between 0 and 1
                //
                F_franjo = 1-C;
     
                // sum all values of F to calculate averages
                Pop[colony_i].mean_F += Pop[colony_i].MyAnts[ant_i].F;
                mean_F_franjo += F_franjo;

                sumsquares_F += Pop[colony_i].MyAnts[ant_i].F *
                    Pop[colony_i].MyAnts[ant_i].F;

                sumsquares_F_franjo += F_franjo * F_franjo;

                // count this ant as an active one (as it has worked
                // at least once)
                ++n_ants_active;
            } // end if workperiods > 0
        } // end for ant_i

        // calculate average switch rate by dividing by the number of active ants
        Pop[colony_i].mean_switches /= n_ants_active;
        
        // calculate mean workperiods by dividing by the total number of ants
        Pop[colony_i].mean_workperiods /= Pop[colony_i].MyAnts.size();

        // calculate variances
        Pop[colony_i].var_switches = sumsquares_switches / n_ants_active - 
            Pop[colony_i].mean_switches * Pop[colony_i].mean_switches;

        Pop[colony_i].var_workperiods = sumsquares_workperiods / Pop[colony_i].MyAnts.size() -
            Pop[colony_i].mean_workperiods * Pop[colony_i].mean_workperiods;

        Pop[colony_i].mean_F /= n_ants_active;
        Pop[colony_i].mean_F_franjo /= n_ants_active;

        Pop[colony_i].var_F = sumsquares_F / n_ants_active - 
            Pop[colony_i].mean_F * Pop[colony_i].mean_F;

        Pop[colony_i].var_F_franjo = sumsquares_F_franjo / n_ants_active - 
            Pop[colony_i].mean_F_franjo * Pop[colony_i].mean_F_franjo;

    } // end for Colonies

} // end of Calc_F()
//==============================================================================================================================================

// Calculates the fitness of a colony from the population size and the number of ants working on each task
// Only observes the later part of the timesteps in a generation, 
// beyond t > tau
// as there is an initialisation effect
void CalcFitness(Population & Pop, Params & Par)
{
    sum_Fit = 0;

    double total;

    bool ant_is_idle;

    for (unsigned int colony_i = 0; colony_i < Pop.size(); ++colony_i)
    {
        Pop[colony_i].idle = 0;

        total = 0;

        // calculate total work periods
        for (int task_i = 0; task_i < Par.tasks; ++task_i)
        {
            total += Pop[colony_i].last_half_acts[task_i];
        }
        
        if (total == 0)
        {
            Pop[colony_i].fitness =0;
        }
        else if (total > 0)
        {
            Pop[colony_i].fitness =  total;

            for (int task_i = 0; task_i < Par.tasks; ++task_i)
            {
                // now multiply by weighted fitness
                Pop[colony_i].fitness *= 
                    pow((Pop[colony_i].last_half_acts[task_i]/total), 
                            Par.fitness_weights[task_i]); 
            }
        }

        // calculate number of idle ants that have never worked
        for (unsigned int ant_i = 0; 
                ant_i < Pop[colony_i].MyAnts.size(); ++ant_i)
        {
            ant_is_idle = true;        

            for (int task_i = 0; task_i < Par.tasks; ++task_i)
            {
                if (Pop[colony_i].MyAnts[ant_i].countacts[task_i] > 0)
                {
                    ant_is_idle = false;
                    break;
                }
            }

            if (ant_is_idle)
            {
                ++Pop[colony_i].idle;
            }
        }
    }

    // find colony with minimum fitness value in population
    int min_fit = Pop[0].ID;

    for (unsigned int col = 1; col < Pop.size(); ++col)
    {
        if (Pop[col].fitness < Pop[min_fit].fitness)
        {
            min_fit = col;
        }
    }

    // calculate difference in fitness between every colony and the lowest colony
    for (unsigned int i = 0; i < Pop.size(); ++i)
    {
        Pop[i].diff_fit = Pop[i].fitness - Pop[min_fit].fitness;
        sum_Fit += Pop[i].diff_fit;
    }

    // make cumulative fitness distribution
    for (unsigned int i = 0; i < Pop.size(); ++i)
    {
        Pop[i].rel_fit = Pop[i].diff_fit / sum_Fit;

        if (i==0)
        {
            Pop[i].cum_fit = Pop[i].rel_fit;
        }
        else
        {
            Pop[i].cum_fit = Pop[i-1].cum_fit + Pop[i].rel_fit;
        }
    }

} // end of CalcFitness()
//-------------------------------------------------------------------------------

// Draw a parent from the cumulative distribution
// Randomly selects genotype samples from all colonies based on their fitness (fittest more likely to be chosen)
int drawParent(int nCol, Population & Pop)
{
    const double draw = gsl_rng_uniform(rng_global); 
    int cmin=-1, cmax=nCol-1;

    // binary search of cumulative fitness value
    while (cmax - cmin != 1)
    {
        int cmid = (cmax+cmin)/2;
        if (draw < Pop[cmid].cum_fit)
        {
            cmax = cmid;
        }
        else
        {
            cmin = cmid;
        }
    }
    return(cmax);
}// end of drawParent()



//Creates the sexual individuals from a colony.
//Repeats drawParent until all required samples are collected.
void MakeSexuals(Population & Pop, Params & Par)
{
    mySexuals.resize(2 * Par.Col); // number of sexuals needed
    parentCol.resize(mySexuals.size());

    for (unsigned int ind = 0; ind < mySexuals.size(); ++ind)
	{
        //initialize sexuals
        mySexuals[ind].threshold.resize(Par.tasks);
        mySexuals[ind].countacts.clear();
        mySexuals[ind].last_act = Par.tasks;
        mySexuals[ind].curr_act = Par.tasks;

        mySexuals[ind].switches = 0;
        mySexuals[ind].F = 0;
        mySexuals[ind].mated = false;
        mySexuals[ind].mated = false;
        
        
        // draw a parent colony for each sexual
        parentCol[ind] = drawParent(Pop.size(), Pop);

        Inherit(mySexuals[ind], 
                Pop[parentCol[ind]].queen, 
                Pop[parentCol[ind]].male, 
                Par);

	}// end for ind
}

// Create the workers in a new colony.
// Randomly selects a mother and father from the sexual individuals of all colonies and combines their genotypes
void MakeColonies(Population &Pop, Params &Par, int generation)
{
    int mother, father;

    // file to write founders from the last generation
    // to, if simulation may be continued at a later time
    static ofstream lastgen;

    for (unsigned int col = 0; col < Pop.size(); ++col)
	{
        do
        {
            mother = gsl_rng_uniform_int(rng_global, mySexuals.size());
            father = gsl_rng_uniform_int(rng_global, mySexuals.size());
        }
        while
		(mother == father 
         || mySexuals[mother].mated==true 
         || mySexuals[father].mated==true);


        mySexuals[mother].mated = true;
        mySexuals[father].mated = true;

        Pop[col].queen = mySexuals[mother]; 
        Pop[col].male = mySexuals[father]; 

	} // end for Colonies


    // if last generation, write out the founders to a file 
    if (generation - simstart_generation == Par.maxgen - 1) 
    {
        lastgen.open("lastgen.txt");
        lastgen << simpart + 1 << endl;
        lastgen << generation + 1 << endl;

        for (unsigned int col = 0; col < Pop.size(); ++col)
        {
            for (int task = 0; task < Par.tasks; ++task)
            {
                lastgen << Pop[col].male.threshold[task] << endl;	
            }
            for (int task = 0; task < Par.tasks; ++task)
            {
                lastgen << Pop[col].queen.threshold[task] << endl;	
            }
        }
    }
} // end MakeColonies()

// give names to each of the datafiles
void NameDataFiles(
        string & data1, 
        string &data2, 
        string &data3, 
        string &data4, 
        string &dataant)
{
    stringstream tmp;
    tmp << "data_work_alloc_" << simpart << ".txt";
    data1 = tmp.str();

    stringstream tmp2;
    tmp2 << "threshold_distribution_" << simpart << ".txt";
    data2 = tmp2.str();

    stringstream tmp3;
    tmp3 << "f_dist_" << simpart << ".txt";
    data3 = tmp3.str();

    stringstream tmp4;
    tmp4 << "branch.txt";
    data4 = tmp4.str();

    stringstream tmp5;
    tmp5 << "ant_beh_" << simpart << ".txt";
    dataant = tmp5.str();
}

// see whether we need to start the simulation from a
// previous instantiation
void Start_Simulation_From_Previous(Params & Par, Population & Pop) 
{
    // if the lastgen.txt file is present
    // this simulation has already got a start and we need
    // to continue it 
    if (FileExists("lastgen.txt"))
    {
        ifstream inp("lastgen.txt");

        // use this to initiate the queens, kings, etc
        StartFromLast(inp, Par, Pop);
    }
	else 
    {
        simpart = 1;
        simstart_generation = 0; 
    }
}

void Update_Col_Data(
        int step,  // current timestep
        Population & Pop, 
        Params & Par)
{
	if (step >= Par.tau) 
    {
        for (unsigned int col = 0; col < Pop.size(); ++col)
        {
            for (int task = 0; task < Par.tasks; ++task)
            {
                Pop[col].last_half_acts[task] += 
                    Pop[col].workfor[task]/Par.alfa[task];
                
                Pop[col].mean_work_alloc[task] += 
                    Pop[col].workfor[task] / Par.alfa[task];
            }	
        }
    }
}

void Write_Headers(
        ofstream &header_file,
        ofstream &data_1gen,
        Params &Par)
{
    header_file << "Gen" << "\t" 
    << "Col"  << "\t";

    data_1gen << "Time" << ";" 
		<< "Col" << ";";

    for (int task_i = 0; task_i < Par.tasks; ++task_i)
    {
        header_file << "NumActs" << (task_i + 1) << "\t"; 
        data_1gen << "Stim" << (task_i + 1) << "\t"; 
    }
    
    for (int task_i = 0; task_i < Par.tasks; ++task_i)
    {
        header_file <<  "Workalloc" << (task_i + 1) << "\t";
        data_1gen << "Workers" << (task_i + 1) << "\t"; 
    }

    header_file << "Idle"<< "\t" 
    << "Inactive" << "\t"
    <<"Fitness" << "\t";
    
    for (int task_i = 0; task_i < Par.tasks; ++task_i)
    {
        header_file <<  "Endstim" << (task_i + 1) << "\t";
    }

    header_file << "mean_switches" << "\t"
    << "mean_workperiods" << "\t";

    for (int task_i = 0; task_i < Par.tasks; ++task_i)
    {
        header_file <<  "last_half_acts" << (task_i + 1) << "\t";
    }

    header_file << endl;
    
	data_1gen << "Fitness" << ";" 
		<< "Mean_F"<< ";" 
		<< "Mean_F_franjo" << endl;
}

//==========================================================================================================
//write out data for graphs 
//==========================================================================================================
void Write_Col_Data(
        ofstream & mydata, 
        Params & Par, 
        Population & Pop, 
        int gen, 
        int colony)
{
	mydata << gen << "\t" << colony << "\t";

	for (int task = 0; task < Par.tasks; ++task)
    {
		    mydata << Pop[colony].numacts[task] << "\t";
    }
	for (int task = 0; task < Par.tasks; ++task)
    {
		    mydata << Pop[colony].mean_work_alloc[task] << "\t"; 				
    }

	mydata << Pop[colony].idle 
            << "\t" << Pop[colony].inactive 
            << "\t" << Pop[colony].fitness;

	for (int task = 0; task < Par.tasks; ++task)
    {
            mydata << "\t" << Pop[colony].stim[task]; 
    }

    mydata << "\t" << Pop[colony].mean_switches 
            << "\t" << Pop[colony].mean_workperiods;

	for (int task = 0; task < Par.tasks; ++task)
    {
            mydata << "\t" << Pop[colony].last_half_acts[task]; 
    }

    mydata << endl;  
}
//------------------------------------------------------------------------------------------------------

void Write_Thresholds_Spec(
        ofstream & data_thresh,
        ofstream & data_f,
        Params & Par,
        Population & Pop, 
        int gen,
        int colony)
{
	 // output only thresholds of foundresses and mean specialization 
	data_thresh << gen; 

    double p_i[Par.tasks];

    double total_acts = 0;

    // calculate denominator for the specialization measure
    double denomin = 0;

	for (int task = 0; task < Par.tasks; ++task)
    {
        data_thresh << ";" << Pop[colony].male.threshold[task];

        total_acts += Pop[colony].numacts[task];
        p_i[task] = Pop[colony].numacts[task];
    }

	for (int task = 0; task < Par.tasks; ++task)
    {
        data_thresh << ";" << Pop[colony].queen.threshold[task];

        p_i[task] /= total_acts;

        denomin += p_i[task]*p_i[task];
    }

    data_thresh << endl;

	//cout << "denominator :" << denomin << endl;

	data_f << gen <<";" 
            << Pop[colony].mean_F <<";"
            << Pop[colony].mean_F_franjo/denomin << ";" 
	    	<< Pop[colony].mean_switches << ";" 
            << Pop[colony].mean_workperiods << ";"
	    	<< Pop[colony].var_F << ";" 
            << Pop[colony].var_F_franjo << ";"
	    	<< Pop[colony].var_switches << ";" 
            << Pop[colony].var_workperiods << endl;
} 

// keep informed of whether branching occurred
void Write_Branching(ofstream & afile, Params & Par, int yesno)
{
	afile << Par.mutp << ";" 
		<< Par.mutstd << ";" 
		<< Par.timecost << ";" 
		<< yesno << endl;	
}	

// check whether there is specailization
bool Check_Spec(Population & Pop)
{
	bool mybool=false;
	int count_cols=0; // count colonies with less than 0.75 average specialization 
	for (unsigned int col=0; col<Pop.size(); col++)
		{
	        	if(Pop[col].mean_F<0.75) 
				count_cols +=1;	

		}

        if (double(count_cols)/Pop.size() < 0.4) mybool=true; // if less than 40% of colonies have low specialization
	
	return mybool; 
} // end of Check_Spec

//=====================================================================================================
// to use in case you want the simulation to stop if specialization is present 
void StopIfSpec(
        Population & Pop, 
        int generation, 
        Params & Par, 
        ofstream & file1, 
        ofstream &file2, 
        ofstream & file3, 
        ofstream &file4, 
        int eq_steps, 
        string filename)
{ 
    if (Check_Spec(Pop))
    {
        file4.open(filename.c_str());
        int branch = 1;

        for (unsigned int col = 0; col < Pop.size(); ++col)
        {
            for (int task = 0; task < Par.tasks; ++task)
            {
                Pop[col].mean_work_alloc[task] /= eq_steps;
            }
            
            Write_Col_Data(file1, Par, Pop, generation, col);
          
            Write_Thresholds_Spec(file2, file3, Par, Pop, generation, col);
        
        }
    
        Write_Branching(file4, Par, branch);	
        exit (1);
    }
    else if (generation == simstart_generation + Par.maxgen-1)
    {
        file4.open(filename.c_str());
    
        int branch = 0;	
        for (unsigned int col = 0; col < Pop.size(); col++)
        {
            for (int task=0; task<Par.tasks; task++)
            {
                Pop[col].mean_work_alloc[task]/=eq_steps;
            }
                
            Write_Col_Data(file1, Par, Pop, generation, col);
              
                 // output only thresholds of foundresses and mean specialization 
                    
            Write_Thresholds_Spec(file2, file3, Par, Pop, generation, col);
            
        }

        Write_Branching(file4, Par, branch);	
    }

} // end StopIfSpec

int main(int argc, char* argv[])
{
	Params myPars;

	ifstream inp("params.txt");
	myPars.InitParams(inp);
	ShowParams(myPars);
	
    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rng_global = gsl_rng_alloc(T);
    gsl_rng_set(rng_global, myPars.seed);

    // also set the seed to be used by random_shuffle
    srand(myPars.seed);

    // initialize the metapopulation
	Population MyColonies;

    // initiliaze the founders
	InitFounders(MyColonies, myPars);
	
	
	// code that checks which run this is
	// is lastgen.txt present y/n
    Start_Simulation_From_Previous(myPars, MyColonies);

    // the datafiles
	string datafile1, 
           threshold_dist_file_name, 
           specialization_dist_file_name, 
           datafile4, 
           dataants;

    // name the files
	NameDataFiles(datafile1, threshold_dist_file_name, specialization_dist_file_name, datafile4, dataants);

	static ofstream out; 
	static ofstream threshold_dist_output_file; // distribution of thresholds across the population
	static ofstream one_generation_output_file;
	static ofstream specialization_dist_output_file;
    static ofstream branching;
    static ofstream header1;
	static ofstream out_ants;
	
    out.open(datafile1.c_str());
    header1.open("header.txt");
    one_generation_output_file.open("data_1gen.txt");

    // write the headers to various data files
    Write_Headers(header1, one_generation_output_file, myPars);

    // open the remaining data files
	threshold_dist_output_file.open(threshold_dist_file_name.c_str()); // threshold distribution file
	specialization_dist_output_file.open(specialization_dist_file_name.c_str());     // 
	    
	out_ants.open(dataants.c_str()); 

    // evolutionary time
	for (int g = simstart_generation; 
            g < simstart_generation + myPars.maxgen; ++g)
    {
        // initialize all colonies in this generation
        Init(MyColonies, myPars);
        
        double equil_steps=0;

	    myPars.gensdone = g;

        // Updates ant behaviour for each timestep (ecological timescale)
        for (int k = 0; k < myPars.maxtime; k++)
        {
            // update all the ants in each colony
            UpdateAnts(MyColonies, myPars);
			myPars.stepsdone = k;

            // update all stimulus levels of each colony
            UpdateStim(MyColonies, myPars);

            // calculate specialization values
            Calc_F(MyColonies, myPars); 

            // update all the fitness data etc
            Update_Col_Data(k, MyColonies, myPars);

            ++equil_steps;
           
            // last timestep before reproduction
            if (k == myPars.maxtime-1)
            {
                // calculate fitness of the colonies
                CalcFitness(MyColonies, myPars);

                // output the data (only at the start or every 100th
                // generation
                if ((g <= 100 || g % 100==0))
                {
                    for (unsigned int col = 0; 
                            col < MyColonies.size(); ++col)
                    {
                        for (int task=0; task<myPars.tasks; ++task)
                        {
                            MyColonies[col].mean_work_alloc[task] /= 
                                equil_steps;
                        }
                           
                        Write_Col_Data(out, myPars, MyColonies, g, col);
                          
                         // output only thresholds of foundresses and mean specialization 
                        Write_Thresholds_Spec(threshold_dist_output_file, 
                                specialization_dist_output_file, 
                                myPars, 
                                MyColonies, 
                                g, 
                                col);
                         
                    } // end of for colonies
                } // end if generations are right
            } // end if k=maxtime

            if (g == simstart_generation + myPars.maxgen - 1) 
            {
                for (unsigned int col = 0; col < MyColonies.size(); ++col)
                {
                    double p1 = (double) MyColonies[col].numacts[0] / 
                        (MyColonies[col].numacts[0] + 
                         MyColonies[col].numacts[1]);
            
                    double p2 = (double)MyColonies[col].numacts[1] / 
                        (MyColonies[col].numacts[0] + 
                         MyColonies[col].numacts[1]);
            
                    double denomin = p1*p1 + p2*p2;

                    one_generation_output_file << k << ";" << col << ";"; 

                    for (int task=0; task<myPars.tasks; task++)
                    {
                         one_generation_output_file <<MyColonies[col].stim[task] << ";"; 
                    }

                    for (int task=0; task<myPars.tasks; ++task)
                    {
                         one_generation_output_file << MyColonies[col].workfor[task]/
                             myPars.alfa[task] << ";";
                    }

                    one_generation_output_file << MyColonies[col].fitness << ";" << 
                            MyColonies[col].mean_F << ";" << 
                            MyColonies[col].mean_F_franjo/denomin << endl; 

                    if (k == myPars.maxtime -1)
                    {	
                        for (unsigned int ant=0; 
                                ant < MyColonies[col].MyAnts.size(); ++ant)
                        {
                            out_ants << col << ";" 
                                << ant << ";" 
                                << MyColonies[col].MyAnts[ant].countacts[0] << ";" 
                                << MyColonies[col].MyAnts[ant].countacts[1] << ";"
                                << MyColonies[col].MyAnts[ant].switches << ";"
                                << MyColonies[col].MyAnts[ant].workperiods << endl;  
                        }
                    }
                }// endfor (unsigned int col
            } // end if if (g == simstart_generation + myPars.maxgen - 1)   
        } // end for (int k = 0; k < myPars.maxtime

        MakeSexuals(MyColonies, myPars);
        MakeColonies(MyColonies, myPars, g);

    } // end for (int g = simstart_generation;
} // end of main()
