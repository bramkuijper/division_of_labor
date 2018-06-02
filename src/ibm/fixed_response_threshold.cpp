
// with stimulus update per ant, or per timestep (see define)
// 02/05/2018 With stochsine function (ctrl+F: //stochsine).
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

//#define DEBUG
#define SIMULTANEOUS_UPDATE

//---------------------------------------------------------------------------

using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *rng_global; // gnu scientific rng 

// Defining parameters
struct Params
{
    int N; //number of workers
    int Col; // number of colonies
    int maxtime; // time steps
    double p;  // quitting probability
    int tasks;
    double mutp;//mutation probability
    int maxgen;
    double beta_fit, gamma_fit;
    int seed;
	
	//stochsine
	double A; //Deterministic factor
	double B; //Stochastic factor
	int genspercycle; //Generations per environmental cycle
	int randommax; //Maximum value of positive random number
	int gensdone; //Generations completed
	int stepsdone; //Timesteps completed in current generation
	double delta;

    int no_task; // number indicating that the current task
                    // is no task

    // vectors containing stimulus parameters
    vector<double> meanT; // mean thresholds
    vector<double> alfa; // stimulus decrease due to work

    // function to read in the parameters from stream
    istream & InitParams(istream & inp);
};


// Defining an individual ant, and the values it contains
struct Ant
{
    // genome
    vector < double > threshold; // thresholds of this ant
    vector < int > countacts;   // counter of acts done by this ant
    int last_act;
    int curr_act;
    int switches; // number of transitions to a different task
    int workperiods; // number of working periods 
    double F; // specialization value
    bool mated;
    double F_franjo;
};

// Vectors of workers and sexuals
typedef vector < Ant > Workers;
typedef vector <Ant> Sexuals;

// Defining the colony of ants, and the values it contains
struct Colony
{
    Workers MyAnts; // stack of workers
    Ant male, queen; // queen and her male
    int ID; // unique ID of the colony

    vector<double> stim; // stimulus levels
    vector<double> newstim; // stimulus levels new timestep
    vector<double> workfor;  // number acts * eff each time step
    vector < int > numacts; // number of acts performed per task
    double idle; // proportion idle workers
    vector <double>last_half_acts; //number of acts performed in the last half of simulation
    double fitness; 
    double diff_fit;// fitness difference to minimal fitness
    double rel_fit; // fitness relative to whole population
    double cum_fit; //cumulative fitness

    vector<double> HighF, LowF, CategF;
    vector<double> HighT, LowT, CategT1, CategT2;
    vector<double>mean_work_alloc;
    double mean_F;
    int num_offspring;
    double mean_F_franjo;
};

// Vector of colonies as the population
typedef vector < Colony > Population;
vector <int> parentCol;
double sum_Fit;
Sexuals mySexuals;

// A For loop/stream that reads in the parameter file, repeating for the number of tasks.
istream & Params::InitParams(istream & in)
{
    // set the number of tasks
    tasks = 2;    

    // read in population size
    in >> N >> // number of workers
        Col >> // number of colonies
        maxtime;  // number of timesteps work is being done
                    // during each evolutionary generation

    double tmp;

    // read in the initial mean thresholds for each stimulus
    for (int i=0; i<tasks; i++) 
    {
        in >> tmp; 
        meanT.push_back(tmp);
    }

    // read in the stimulus decrease parameters
    for (int i=0; i<tasks; i++)
    {
        in >> tmp; 
        alfa.push_back(tmp);
    }

    in >> p >> // quitting probability
        mutp >> // mutation probability
        maxgen >> // maximum number of generations
        beta_fit >> // fitness weights
        gamma_fit >> // fitness weights
		
		//stochsine			 
		A >> //Deterministic factor
		B >> //Stochastic factor
		genspercycle >> //Generations per environmental cycle
		randommax >> //Maximum value of positive random number
		
		seed;

    // set a number which indicates that ants are currently 
    // working on no task at all
    no_task = 7;

    return in;
}
//----------------------------------------------------------------------------------

// Diagnostic function: prints the results of the stream to verify that it performed correctly.
void ShowParams(Params & Par)
{
     cout <<"Workers " << Par.N << endl;
     cout << "Colonies " << Par.Col << endl;
     cout << "Timesteps " << Par.maxtime << endl;
     
     for (int task=0; task<Par.tasks; task++)
     {
        cout << "Initial T" <<task <<"\t" << Par.meanT[task] << endl;
     }
     
     for (int task=0; task<Par.tasks; task++)
     {
        cout << "effic " << task << "\t" << Par.alfa[task] << endl;
     }
     
     cout << "prob quit" << Par.p << endl;
     cout << "Task number " << Par.tasks << endl;
     cout << "mut prob " << Par.mutp << endl;
     cout << "Max gen " <<  Par.maxgen << endl;
     cout << "Exp task 1 " <<  Par.beta_fit << endl;
     cout << "Exp task 2 " <<  Par.gamma_fit << endl;

	   //stochsine
	   cout << "Deterministic Factor A " << Par.A << endl;
	   cout << "Stochastic Factor B " << Par.B << endl;
	   cout << "Generations per Cycle " << Par.genspercycle << endl;
	   cout << "Maximum Random Number " << Par.randommax << endl;

	   cout << "Seed " << Par.seed << endl;
}
//-----------------------------------------------------------------------------

// Initialises the founding generation by mating a male and a queen, with attributes randomly selected from a normal distribution
void InitFounders(Population &Pop, Params &Par)
{
#ifdef DEBUG  
    assert(Par.tasks==2);
    cout <<Par.Col << endl;
#endif
    Pop.resize(Par.Col);

    // initialize thresholds of males and queens
    for (unsigned int i = 0; i < Pop.size(); ++i)
    {
        Pop[i].male.threshold.resize(Par.tasks);
        Pop[i].queen.threshold.resize(Par.tasks);

        // give males and queens random thresholds
        for (int task=0; task<Par.tasks; ++task)
        {
            Pop[i].male.threshold[task] = Par.meanT[task] 
                + gsl_ran_gaussian(rng_global,1);

            Pop[i].queen.threshold[task] = Par.meanT[task] 
                + gsl_ran_gaussian(rng_global,1);
        }

        Pop[i].male.mated =true;
        Pop[i].queen.mated = true;
    }
}
//----------------------------------------------------------------------------------------------------------------------

// Defines inheritance, producing worker threshold genotypes from their parents, including mutation
void Inherit(Ant &Daughter, Ant &Mom, Ant &Dad, Params &Par)
{
    for (int task = 0; task < Par.tasks; ++task)
    {
        // draw a random number designating the parent who is the one
        // who inherits the trait
        int from_which_parent = gsl_rng_uniform_int(rng_global, 2);

        // inherit from mom
        if (from_which_parent == 0)
        {
            // mutate 
            if (gsl_rng_uniform(rng_global) < Par.mutp)
            {
                Daughter.threshold[task] = Mom.threshold[task] 
                    + gsl_ran_gaussian(rng_global,1);
            }
            else 
            {
                Daughter.threshold[task] = Mom.threshold[task];
            }
        }
        else // alternatively, inherit from dad
        {
            if (gsl_rng_uniform(rng_global) < Par.mutp)
            {
                Daughter.threshold[task] = Dad.threshold[task] 
                    + gsl_ran_gaussian(rng_global,1);
            }
            else
            {
                Daughter.threshold[task] = Dad.threshold[task];
            }
        }
        
        // perform boundary checking, as thresholds cannot be
        // negative
        if (Daughter.threshold[task] < 0)
        {
            Daughter.threshold[task] = 0.0;
        }
    }
}
//----------------------------------------------------------------------------------------------------------

// Initialises the ant workers, as generated by the Inherit function
// Prints their threshold variables, and puts them to work on a task depending on their inherited thresholds
void InitAnts(Ant & myAnt, Params & Par, Colony & myCol)
{
    myAnt.threshold.resize(Par.tasks);

    if (Par.maxgen > 1)
    {
        Inherit(myAnt, myCol.queen, myCol.male, Par);
    }
    else 
    {
        for (int task=0; task<Par.tasks; ++task)
        {
            myAnt.threshold[task] = Par.meanT[task];
        }
    }

    myAnt.countacts.resize(Par.tasks);

    for (int task = 0; task < Par.tasks; ++task)
    {
        myAnt.countacts[task]=0;
    }
    
    myAnt.last_act = Par.no_task; // set at a value of no task
    myAnt.curr_act = Par.tasks; // 0 = task 1, 1 = task 2, etc, n = idle
    myAnt.switches = 0; // set counter of switches to 0
    myAnt.workperiods=0;
    myAnt.F = 10;
    myAnt.F_franjo = 10;
    myAnt.mated = false;
}
//------------------------------------------------------------------------------------

// Initilializes all colonies at the start of each evolutionary generation
void Init(Population & Pop, Params & Par)
{
#ifdef DEBUG
    cout << Pop.size() << endl;
    cout << Par.N << endl;
    assert(Pop.size()==Par.Col);
#endif
      
    for (unsigned int i = 0; i< Pop.size(); i++)
    {
        Pop[i].MyAnts.resize(Par.N);
        Pop[i].ID = i;
        Pop[i].fitness = 0;
        Pop[i].rel_fit = 0;
        Pop[i].cum_fit = 0;
        Pop[i].diff_fit = 0;
        Pop[i].idle = 0;
        Pop[i].mean_F = 10; // specialization value
        Pop[i].mean_F_franjo =10;

        // erase any existing workforce numbers
        if (!Pop[i].workfor.empty())
        {
            Pop[i].workfor.erase(
                    Pop[i].workfor.begin(),
                    Pop[i].workfor.end());
        }

        assert(Pop[i].workfor.empty());

        // erase any existing fitness data
        if (!Pop[i].last_half_acts.empty()) 
        {
            Pop[i].last_half_acts.erase(
                    Pop[i].last_half_acts.begin(),
                    Pop[i].last_half_acts.end());
        }

        assert(Pop[i].last_half_acts.empty()); 

        if (!Pop[i].numacts.empty()) 
        {
            Pop[i].numacts.erase(
                    Pop[i].numacts.begin(),
                    Pop[i].numacts.end());
        }

        assert(Pop[i].numacts.empty()); 

        if (!Pop[i].stim.empty()) 
        {
            Pop[i].stim.erase(
                    Pop[i].stim.begin(),
                    Pop[i].stim.end());
        }

        assert(Pop[i].stim.empty()); 

        if (!Pop[i].newstim.empty())
        {
            Pop[i].newstim.erase(
                    Pop[i].newstim.begin(),
                    Pop[i].newstim.end());
        }

        assert(Pop[i].newstim.empty());

        if (!Pop[i].mean_work_alloc.empty())
        {
            Pop[i].mean_work_alloc.erase(
                    Pop[i].mean_work_alloc.begin(),
                    Pop[i].mean_work_alloc.end());
        }

        assert( Pop[i].mean_work_alloc.empty());

        for (int job = 0; job < Par.tasks; ++job)
        {
            Pop[i].workfor.push_back(0);
            Pop[i].last_half_acts.push_back(0);
            Pop[i].numacts.push_back(0);
            Pop[i].stim.push_back(0);
            Pop[i].newstim.push_back(0);
            Pop[i].mean_work_alloc.push_back(0);
        }
    
        for (unsigned int j = 0; j < Pop[i].MyAnts.size(); ++j)
        {
            InitAnts(Pop[i].MyAnts[j], Par, Pop[i]);
        }
    } // end of for Pop
} // end of Init()
//-------------------------------------------------------------------------------

// Updates the stimulus per ant for each timestep.
// Stimulus for a task decreased based of how much work put towards it.
void UpdateStimPerAnt(Params & Par, Colony & anyCol, Ant & anyAnt, int task)
{
#ifdef DEBUG
    cout << Par.alfa[task] << endl;
    cout << Par.N << endl;
    cout << anyCol.workfor[task] << endl;
    cout <<anyCol.numacts[task]<< endl;
    cout << anyCol.stim[task] << endl;
#endif
    anyCol.stim[task] -= (Par.alfa[task]/Par.N); 

    if (anyCol.stim[task] < 0)
    {
        anyCol.stim[task] = 0;
    }
}
//------------------------------------------------------------------------------

// Checks whether an ant will cease task performance and becomes idle
// Based on random draws of an ant's innate quitting probability, unrelated to any other varaiable
void QuitTask(Colony & anyCol, Ant & anyAnt, int job, Params & Par)
{
#ifdef DEBUG
cout << "Quitting tasks" << endl;
cout << "chance to quit: " << Par.p << endl;
#endif

    // ant quits
    if (gsl_rng_uniform(rng_global) < Par.p)
    {
        anyAnt.curr_act = Par.tasks;
    }
#ifndef SIMULTANEOUS_UPDATE  
    else
    {
        UpdateStimPerAnt(Par, anyCol, anyAnt, job);
    }
#endif
}
//------------------------------------------------------------------------------

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
    else if (s > 0.00000)
    {
        RP = 1.0; // if threshold is zero, then probability of acting is 1.  
    }
    else
    {
        RP=0.0;
    }
    
    return(RP);
}
//-------------------------------------------------------------------------------

// Causes an ant to perform one of the tasks, using responce probability
// Involves the stimulus per ant and threshold, and defines when an ant might switch tasks
void TaskChoice(Params & Par, Colony & anyCol, Ant & anyAnt)
{ 
    // prob to meet one of two tasks is random
    int job =  gsl_rng_uniform_int(rng_global, Par.tasks);  
       
#ifdef DEBUG 
   cout << "Choosing tasks" << endl; 
    assert(Par.tasks==2);
    assert(anyCol.stim[job]>=0);
    assert(anyAnt.threshold[job]>=0);
#endif

    // ant chooses to perform the randomly chosen job
    if (gsl_rng_uniform(rng_global) < 
            RPfunction(anyAnt.threshold[job],anyCol.stim[job])) 
    {
        anyAnt.curr_act = job; 
        anyAnt.workperiods += 1; 

#ifndef SIMULTANEOUS_UPDATE 
         UpdateStimPerAnt(Par, anyCol, anyAnt, job);
#endif

    } 
    else //randomly encountered task not chosen, ant stays idle
    {
        anyAnt.curr_act=2;
    }
} // end of TaskChoice()
//-------------------------------------------------------------------------------------------------

// Reassesses whether an ant should be working or idle
// Updates work done towards a task using QuitTask and TaskChoice on all ants per colony
void UpdateAnts(Population & Pop, Params & Par)
{
    for (unsigned int colony_i = 0; colony_i < Pop.size(); ++colony_i)
    {
        for (int task = 0; task < Par.tasks; ++task)
        {
            // reset work for tasks
            Pop[colony_i].workfor[task] = 0; 
        }

        // let active ants potentially quit
        // let idle ants potentially find work
        for (unsigned int ant_i = 0; 
                ant_i < Pop[colony_i].MyAnts.size(); ++ant_i)  
        {
            //if ant currently active 
            if (Pop[colony_i].MyAnts[ant_i].curr_act < Par.tasks)
            {
                //record last act
                Pop[colony_i].MyAnts[ant_i].last_act = 
                    Pop[colony_i].MyAnts[ant_i].curr_act; 
                
                // ant wants to quit?
                QuitTask(Pop[colony_i], 
                        Pop[colony_i].MyAnts[ant_i],
                        Pop[colony_i].MyAnts[ant_i].curr_act,
                        Par); 
            }
            else // ant currently not active
            {
                //if inactive, choose a task 
                TaskChoice(Par, Pop[colony_i], Pop[colony_i].MyAnts[ant_i]); 
            }

            //update number of switches after choosing tasks
            //
            // ant is currently active
            if (Pop[colony_i].MyAnts[ant_i].curr_act < Par.tasks)
            {
                // get current act 
                int current_act_id = Pop[colony_i].MyAnts[ant_i].curr_act;

                // update colony-level counter of number of acts on the current task
                ++Pop[colony_i].numacts[current_act_id];

                // update ant-level counter of number of acts on the current task
                ++Pop[colony_i].MyAnts[ant_i].countacts[current_act_id];

                // update the effective amount work done on the task
                Pop[colony_i].workfor[current_act_id] += Par.alfa[current_act_id]; 
                
                // update the work done for that task 
                if (Pop[colony_i].MyAnts[ant_i].last_act < Par.tasks 
                        && Pop[colony_i].MyAnts[ant_i].last_act != 
                                Pop[colony_i].MyAnts[ant_i].curr_act)
                {
                    ++Pop[colony_i].MyAnts[ant_i].switches;
                }
            }
        } // end for (unsigned int ant_i 
    } // end for (unsigned int colony_i 
}  // end of UpdateAnts()
//------------------------------------------------------------------------------

//stochsine
//Creates value for delta with a stochastic sine wave (Botero et al. 2015)
void Stochsine(Params & Par)
{
Par.delta = (Par.A  * sin((2 *

	//pi
	3.14159265358979323846 *

	//Calculate cumulative timesteps
	((Par.maxtime*Par.gensdone) + Par.stepsdone)

	) / Par.maxtime* Par.genspercycle))
	+ (Par.B *

	//Random number between 0 and randdommax
	gsl_rng_uniform_pos(rng_global) * Par.randommax);
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
        // update the stimulus for each task
        for (int task = 0; task < Par.tasks; ++task)
        {
            // calculate new stimulus level
            Pop[colony_i].newstim[task] = Pop[colony_i].stim[task] 
                + Par.delta;

#ifdef SIMULTANEOUS_UPDATE
            Pop[colony_i].newstim[task] -= (Pop[colony_i].workfor[task]/Par.N);
#endif

            // update the stimulus
            Pop[colony_i].stim[task] = Pop[colony_i].newstim[task];

            if (Pop[colony_i].stim[task] < 0)
            {
                Pop[colony_i].stim[task] = 0;
            }
        }
    }
}

//------------------------------------------------------------------------------

// Calculates specialization
// Involves counting the number of task switches and the number of work periods (acts)
void Calc_F(Population & Pop, Params & Par)
{
    double C;
    double totacts; // count of total acts on a task per ant

    for (unsigned int colony_i = 0; 
            colony_i < Pop.size(); 
            ++colony_i)
    {
        
        
        double sumF=0;
        double sumF_franjo=0;
        double activ=0;

        for (unsigned int ant_i = 0; 
                ant_i < Pop[colony_i].MyAnts.size(); 
                ++ant_i)
        {
            totacts = 0;
            C = 0;
            
            // sum the total amount acts over all tasks
            for (int task = 0; task < Par.tasks; ++task)
            {
                totacts += Pop[colony_i].MyAnts[ant_i].countacts[task];
            }

            // calculate specialization values per ant
            if (totacts > 0) 
            {
                // the probability to switch = switches / acts
                C = double(Pop[colony_i].MyAnts[ant_i].switches) / 
                    Pop[colony_i].MyAnts[ant_i].workperiods;

                // q in Duarte et al 2012 BES eq. (5) is 1 - C

                Pop[colony_i].MyAnts[ant_i].F = 1.0 - 2.0*C;
                Pop[colony_i].MyAnts[ant_i].F_franjo = 1.0 - C;
            }
            
            if (totacts > 0 && 
                    Pop[colony_i].MyAnts[ant_i].curr_act < Par.tasks)
            {
                sumF += Pop[colony_i].MyAnts[ant_i].F;
        
                activ +=1;

                sumF_franjo += Pop[colony_i].MyAnts[ant_i].F_franjo;
            }
        } // end for (unsigned int ant_i = 0; 

        Pop[colony_i].mean_F = sumF/activ;
        Pop[colony_i].mean_F_franjo = sumF_franjo/activ; 
    }
}
//------------------------------------------------------------------------------

// Calculates the fitness of a colony from the population size and the number of ants working on each task
// Only observes the second half of the timesteps in a generation, as there is an initialisation effect
void CalcFitness(Population & Pop, Params & Par)
{
    sum_Fit = 0;
    double total = 0;

    for (unsigned int colony_i = 0; colony_i < Pop.size(); ++colony_i)
    {
        Pop[colony_i].idle = 0;

        total = Pop[colony_i].last_half_acts[0] 
            + Pop[colony_i].last_half_acts[1];

        if (total == 0)
        {
            Pop[colony_i].fitness =0;
        }
        else if (total > 0)
        {
            Pop[colony_i].fitness =  total * 
                pow((Pop[colony_i].last_half_acts[0]/total), Par.beta_fit) * 
                pow((Pop[colony_i].last_half_acts[1]/total), Par.gamma_fit);
        }

        for (unsigned int ant_i = 0; 
                ant_i < Pop[colony_i].MyAnts.size(); ++ant_i)
        {
            if (Pop[colony_i].MyAnts[ant_i].countacts[0] == 0 
                    && Pop[colony_i].MyAnts[ant_i].countacts [1] == 0)
            {
                ++Pop[colony_i].idle;
            }
        }
    }

    int min_fit = Pop[0].ID;

    for (unsigned int col = 1; col < Pop.size(); ++col)
    {
        if (Pop[col].fitness < Pop[min_fit].fitness)
        {
            min_fit = col;
        }
    }

    for (unsigned int i = 0; i < Pop.size(); ++i)
    {
        Pop[i].diff_fit = Pop[i].fitness - Pop[min_fit].fitness;
        sum_Fit += Pop[i].diff_fit;
    }

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

// Puts workers into categories of low and high response thresholds for data analysis
void Categorize(Population & Pop)
{
    for (unsigned int colony_i = 0; colony_i < Pop.size(); ++colony_i)
    {   //creating bins
        Pop[colony_i].LowF[0] = -1;
        Pop[colony_i].HighF[0] = -0.9;
        Pop[colony_i].CategF[0] = 0;
        Pop[colony_i].LowT[0] = 0;
        Pop[colony_i].HighT[0] = 0.5;
        Pop[colony_i].CategT1[0] = 0;
        Pop[colony_i].CategT2[0]=0;

        for (unsigned int index = 1; index < Pop[colony_i].HighF.size(); index++)
        {
            Pop[colony_i].LowF[index] = Pop[colony_i].LowF[index-1]+0.1;
            Pop[colony_i].HighF[index] = Pop[colony_i].HighF[index-1]+0.1;
            Pop[colony_i].CategF[index] = 0;

            Pop[colony_i].LowT[index]= Pop[colony_i].LowT[index-1]+0.5;
            Pop[colony_i].HighT[index] = Pop[colony_i].HighT[index-1]+0.5;
            Pop[colony_i].CategT1[index] = 0;
            Pop[colony_i].CategT2[index] = 0;
        }

        // TODO
        Pop[colony_i].LowF[10] = 0;
        Pop[colony_i].HighF[9] = 0;

        // relative frequencies
        if (Pop[colony_i].MyAnts.size() - Pop[colony_i].idle > 0)
        {
            for (unsigned int worker = 0; 
                    worker < Pop[colony_i].MyAnts.size(); ++worker)
            {
                for (unsigned int index = 0; 
                        index < Pop[colony_i].CategF.size(); ++index)
                {
                    // specialization
                    if (Pop[colony_i].MyAnts[worker].F >= Pop[colony_i].LowF[index]
                            && Pop[colony_i].MyAnts[worker].F < Pop[colony_i].HighF[index])
                    {
                        Pop[colony_i].CategF[index] += 1.0/(
                                (Pop[colony_i].MyAnts.size())-(Pop[colony_i].idle)
                                );
                    }

                    //threshold 1
                    if (Pop[colony_i].MyAnts[worker].threshold[0] >= 
                            Pop[colony_i].LowT[index]
                        && 
                        Pop[colony_i].MyAnts[worker].threshold[0] < 
                            Pop[colony_i].HighT[index])
                    {
                        ++Pop[colony_i].CategT1[index];
                    }

                    // threshold 2
                    if (Pop[colony_i].MyAnts[worker].threshold[1] >= 
                            Pop[colony_i].LowT[index]
                        && 
                        Pop[colony_i].MyAnts[worker].threshold[1] < 
                        Pop[colony_i].HighT[index])
                    {
                        ++Pop[colony_i].CategT2[index];
                    }
                }  // for category
            } // for worker
        } // end if
    }// end of for colony
} // end of CategorizeF()

//------------------------------------------------------------------------------

// Draw a parent from the cumulative distribution
// Randomly selects genotype samples from all colonies based on their fitness (fittest more likely to be chosen)
int drawParent(int nCol, Population & Pop)
{
    const double draw = gsl_rng_uniform(rng_global); 

    int cmin = -1;
    int cmax = nCol - 1;

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
    return cmax;
}// end of drawParent()
//----------------------------------------------------------------------------------------------------

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

        mySexuals[ind].countacts.resize(0);
        mySexuals[ind].last_act = Par.tasks;
        mySexuals[ind].curr_act = Par.tasks;
        mySexuals[ind].switches = Par.tasks;
        mySexuals[ind].F = 0;
        mySexuals[ind].mated = false;

        // draw a parent colony for each sexual
        parentCol[ind] = drawParent(Pop.size(), Pop);

        Inherit(mySexuals[ind], 
                Pop[parentCol[ind]].queen, 
                Pop[parentCol[ind]].male, 
                Par);

    }
}
//-------------------------------------------------------------------------------------------

// Create the workers in a new colony.
// Randomly selects a mother and father from the sexual individuals of all colonies and combines their genotypes
void MakeColonies(Population &Pop, Params &Par)
{
    int mother, father;

    for (unsigned int col = 0; col < Pop.size(); ++col)
    {
        do
        {
            mother = gsl_rng_uniform_int(rng_global, mySexuals.size());
            father = gsl_rng_uniform_int(rng_global, mySexuals.size());
        }
        while (mother == father || 
                mySexuals[mother].mated || 
                mySexuals[father].mated);

        mySexuals[mother].mated = true;
        mySexuals[father].mated = true;

        Pop[col].queen = mySexuals[mother]; 
        Pop[col].male = mySexuals[father]; // copy new males and females

    } // end for Colonies
} // end MakeColonies()
//-----------------------------------------------------------------------------------------------------

//Defines the auxiliary function compare_vector_elements, which sorts vector elements.
//Used because the default C++ sorting function is inflexible.
bool compare_vector_elements(vector<int>anyvector, int size_vector)
{
    bool mybool=false;

    for (int m = 0; m < size_vector-1; ++m) 
    {
        for (int n = m+1; n < size_vector; ++n)
        {
            if (anyvector[m] == anyvector[n])
            {
                mybool=false;
            
                break;
            }
            else
            {
                mybool=true;
            }
        }
    }

    return(mybool);
}
//------------------------------------------------------------------------------------------------------

//The main body of the program.
//Starts by creating the output streams for the results.
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
    
    static ofstream out; 
    static ofstream out2;
    static ofstream out3;
    static ofstream out4;

    if(myPars.Col>1) 
    {
        out.open("data_work_alloc.txt");
        out3.open("stimulus_acts.txt");
        
        out << "Gen" << "\t" << "Col"  << "\t" << "A"  << "\t" << "B"  << "\t" << "NumActs1" << "\t" << "NumActs2" <<
        "\t" << "WorkAlloc1" << "\t" << "WorkAlloc2" <<"\t" << "Idle"<< "\t" << "Fitness" << "\t" << "Mean_F" << "\t" << "Mean_F_franjo" <<endl; 

        out3 << "Gen" << ";" << "Time" << ";" << "Col" << ";" << "Stim1" << ";" << "Workers1" << ";" << "Stim2" << ";" << "Workers2" << ";" << "Fitness" << endl;
    }
    else 
    {
        out4.open("data_1gen.txt");

        out4 << "Time" << ";" << "Col" << ";" << "Stim1" << ";" << "Stim2" << ";" << "Workers1" << ";" << "Workers2" << ";" << "Fitness" << ";" 
             << "Mean_F" << ";" << "Mean_F_franjo" << endl;
    } 

    out2.open("threshold_distribution.txt");
    
    
    //headers for data
    Population MyColonies;
    InitFounders(MyColonies, myPars);

    for (int g = 0; g <= myPars.maxgen; ++g)
    {
        Init(MyColonies, myPars);
        
        double equil_steps=0;

        // Picks ten random colonies as founders and shuffles them:
        vector<int>random_colonies(10);

#ifdef DEBUG
    cout << MyColonies.size() <<endl;
    cout <<myPars.maxtime << endl;
#endif

        // Updates ant behaviour for each timestep (ecological timescale)
        for (int k = 0; k < myPars.maxtime; ++k)
        {
            UpdateAnts(MyColonies, myPars);
			//stochsine
			myPars.gensdone = g;
			myPars.stepsdone = k;
            UpdateStim(MyColonies, myPars);
       
            if (k >= myPars.maxtime/2)
            {
                for (unsigned int col = 0; col < MyColonies.size(); ++col)
                {
                    for (int task=0; task<myPars.tasks; ++task)
                    {
                        MyColonies[col].last_half_acts[task] += 
                            MyColonies[col].workfor[task]/myPars.alfa[task];
                    }
                }	
            }

            // Calculate work allocation for the last fifty timesteps
            if (k >= myPars.maxtime-51) 
            {
                equil_steps +=1;

                for (unsigned int col = 0; col < MyColonies.size(); col++)
                {
                    for (int task=0; task<myPars.tasks; task++)
                    {
                        MyColonies[col].mean_work_alloc[task] += 
                            MyColonies[col].workfor[task]/myPars.alfa[task];
                    }
                }
            }
            
			
			// Calculates fitness, categorises ants 
            // into low/high threshold categories, giving them as an output
            Calc_F(MyColonies, myPars);



            // write the values of the stimulus levels for each colony 
            // to "stimulus_acts.txt"
            //
            // only output stimulus every nth timestep to prevent datafiles becoming
            // massive. If you want to output it every timestep, set
            // k % 1
            if (k % 10 == 0 && myPars.Col > 1)
            {
                // loop through colonies
                for (unsigned int col = 0; col < MyColonies.size(); ++col)
                {
                    // write generation, timestep and colonynumber to file
                    out3 << g << ";" 
                        << k << ";" 
                        << col << ";";

                    // write stimulus levels and worker numbers to file
                    for (int task = 0; task < myPars.tasks; ++task)
                    {
                        out3 << MyColonies[col].stim[task] << ";"
                                << MyColonies[col].workfor[task] << ";";
                    }

                    // calculate colony level fitness and write that to a file
                    out3 << MyColonies[col].fitness << endl;
                }

            } // done writing stimulus and fitness 


            if (k == myPars.maxtime-1)
            {
                CalcFitness(MyColonies, myPars);

                if ((g <= 100 || g % 10 == 0) && myPars.Col > 1)
                {
                    for (unsigned int col = 0; col < MyColonies.size(); ++col)
                    {
                        for (int task = 0; task < myPars.tasks; ++task)
                        {
                            MyColonies[col].mean_work_alloc[task]/=equil_steps;
                        }

                        out << g << "\t" << col << "\t" << myPars.A << "\t" << myPars.B << "\t"; 

                        for (int task = 0; task < myPars.tasks; ++task) 
                        {
                            out << MyColonies[col].numacts[task] << "\t";
                        }

                        for (int task = 0; task < myPars.tasks; ++task)
                        {
                            out << MyColonies[col].mean_work_alloc[task] << "\t"; 
                        }

                        out << "\t" << MyColonies[col].idle << 
                            "\t" << MyColonies[col].fitness << 
                            "\t" << MyColonies[col].mean_F <<
                            "\t" << MyColonies[col].mean_F_franjo << endl;
                 
                        
                        // output only thresholds of foundresses and mean specialization 
                        out2 << g << ";";

                        for (int task=0; task<myPars.tasks; ++task)
                        {
                            out2 << MyColonies[col].male.threshold[task] << ";"; 
                        }

                        for (int task=0; task<myPars.tasks; ++task)
                        {
                            out2 << MyColonies[col].queen.threshold[task] << ";"; 
                        }

                        out2 << MyColonies[col].mean_F <<";" 
                            << MyColonies[col].mean_F_franjo << endl;
                            
                    }
                } // if (g<=100...
            }
             
			// Reduced version of the code that only runs for 1 colony, to test the simulation without high system requirements
            if (myPars.Col == 1)
            {
                for (unsigned int col = 0; col < MyColonies.size(); ++col)
                {
                    out4 << k << ";" << col << ";"; 

                    for (int task=0; task<myPars.tasks; ++task)
                    {
                        out4 <<MyColonies[col].stim[task] << ";"; 
                    }

                    for (int task=0; task<myPars.tasks; ++task)
                    {
                        out4 << MyColonies[col].workfor[task]/myPars.alfa[task] << ";";
                    }

                    out4 << MyColonies[col].fitness << ";" 
                        << MyColonies[col].mean_F << ";" 
                        << MyColonies[col].mean_F_franjo << endl;	     
                    
                    for (unsigned int ind = 0; ind < MyColonies[col].MyAnts.size(); ++ind)
                    {
                        out2 << g << ";";

                        for (int task=0; task<myPars.tasks; ++task)
                        {
                            out2 << MyColonies[col].MyAnts[ind].threshold[task] << ";"; 

                            out2 << MyColonies[col].MyAnts[ind].F << ";" 
                                << MyColonies[col].MyAnts[ind].F_franjo <<endl;
                        }
                    }
                }
            }
        } // for (int k = 0; k < myPars.maxtime; ++k)

        // now make sexuals and new colonies
        MakeSexuals(MyColonies, myPars);
        MakeColonies(MyColonies, myPars);

    } // end for (int g = 0; g < myPars.maxgen; ++g)
} // end main
