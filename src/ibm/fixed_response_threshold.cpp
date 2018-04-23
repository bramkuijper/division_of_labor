
// with stimulus update per ant, or per timestep (see define)

//---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "random.h"
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

// parameter
struct Params
{
    int N; //number of workers
    int Col; // number of colonies
    int maxtime; // time steps
    double p;  // quitting probability
    int tasks;
    double mutp;//mutation probability
    int maxgen;
    double beta_fit, gama_fit;
    int seed;

    int no_task; // number indicating that the current task
                    // is no task

    // vectors containing stimulus parameters
    vector<double> meanT; // mean thresholds
    vector<double> delta; // stimulus increase
    vector<double> alfa; // stimulus decrease due to work

    // vectors containing fitness parameters
    vector<double> beta; // weights of work on task i 
                        // in fitness equation
   
    // function to read in the parameters from stream
    istream & InitParams(istream & inp);
};


// the ant's data
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

// set up vectors of workers and sexuals
typedef vector < Ant > Workers;
typedef vector <Ant> Sexuals;

// set up the colony containing ants
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

    vector<double>HighF, LowF, CategF;
    vector<double> HighT, LowT, CategT1, CategT2;
    vector<double>mean_work_alloc;
    double mean_F;
    int num_offspring;
    double mean_F_franjo;
};

typedef vector < Colony > Population;
vector <int> parentCol;
double sum_Fit;
Sexuals mySexuals;


// read in the parameter file
istream & Params::InitParams(istream & in)
{
    // set the number of tasks
    tasks = 2;    

    // read in population size
    in >> N >>
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

    // read in the baseline increases per stimulus
    for (int i=0; i<tasks; i++) 
    {    
        in >> tmp ; 
        delta.push_back(tmp);
    }

    // read in the stimulus decrease parameters
    for (int i=0; i<tasks; i++)
    {
        in >> tmp; 
        alfa.push_back(tmp);
    }

    // read in the fitness weights of each task
    for (int i=0; i<tasks; i++) 
    {
        in >> tmp;
        beta.push_back(tmp);
    }

    in >> p >> // quitting probability
        mutp >> // mutation probability
        maxgen >> // maximum number of generations
        beta_fit >> // fitness weights
        gama_fit >> // fitness weights
        seed;

    // set a number which indicates that ants are currently 
    // working on no task at all
    no_task = 7;

    return in;
}

//----------------------------------------------------------------------------------

// diagnostic function to see whether parameters
// have been read in correctly
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
        cout << "Delta " <<task <<"\t"  << Par.delta[task] << endl;
     }
     
     for (int task=0; task<Par.tasks; task++)
     {
        cout << "effic " << task << "\t" << Par.alfa[task] << endl;
     }
     
     for (int task=0; task<Par.tasks; task++)
     {
        cout << "Decay " << task << "\t" << Par.beta[task] << endl;
     }
     
     cout << "prob quit" << Par.p << endl;
     cout << "Task number " << Par.tasks << endl;
     cout << "mut prob " << Par.mutp << endl;
     cout << "Max gen " <<  Par.maxgen << endl;
     cout << "Exp task 1 " <<  Par.beta_fit << endl;
     cout << "Exp task 2 " <<  Par.gama_fit << endl;
     cout << "seed " << Par.seed << endl;
}
//-----------------------------------------------------------------------------


// initialize founders 
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

// parent offspring inheritance
void Inherit(Ant &Daughter, Ant &Mom, Ant &Dad, Params &Par)
{
    for (int task = 0; task < Par.tasks; ++task)
    {
        // draw a random number designating the parent who is the one
        // who inherits the trait
        int from_which_parent = gsl_rng_uniform_int(ran_gaussian, 2);

        // inherit from mom
        if (from_which_parent == 0)
        {
            // mutate 
            if (gsl_rng_uniform(ran_gaussian) < Par.mutp)
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
            if (Par.mutp > anynumber)
            {
                Daughter.threshold[task] = Dad.threshold[task] 
                    + gsl_ran_gaussian(rng_global,1);
            }
            else
            {
                Daughter.threshold[task] = Dad.threshold[task];
            }
        }
        
        // set boundaries 
        if (Daughter.threshold[task] < 0)
        {
            Daughter.threshold[task] = 0.0;
        }
    }
}

//----------------------------------------------------------------------------------------------------------


// initialize the ants in the colony
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
    myAnt.curr_act = 2; // 0 = task 1, 1 = task 2, 2 = idle
    myAnt.switches = 0; // set counter of switches to 0
    myAnt.workperiods=0;
    myAnt.F = 10;
    myAnt.F_franjo = 10;
    myAnt.mated = false;
}


//------------------------------------------------------------------------------------

// initilialize all colonies at the start of each 
// evolutionary generation
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
        if (!Pop[i].workfor.empty())
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

// update the stimulus levels per ant
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

// check whether an ant is going to quit
void QuitTask(Colony & anyCol, Ant & anyAnt, int job, Params & Par)
{
#ifdef DEBUG
cout << "Quitting tasks" << endl;
cout << "chance to quit: " << Par.p << endl;
#endif

    // ant quits
    if (gsl_rng_uniform(rng_global) < Par.p)
    {
        anyAnt.curr_act = 2;
    }
#ifndef SIMULTANEOUS_UPDATE  
    else
    {
        UpdateStimPerAnt(Par, anyCol, anyAnt, job);
    }
#endif
}

//------------------------------------------------------------------------------

// response probability
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

// an ant chooses a certain task to perform
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

// update all ants of the colony
void UpdateAnts(Population & Pop, Params & Par)
{
    for (unsigned int colony_i = 0; colony_i < Pop.size(); ++colony_i)
    {
        for (int task = 0; task < Par.tasks; ++task)
        {
            // reset work for tasks
            Pop[colony_i].workfor[task] = 0; 
        }

        for (unsigned int ant_i = 0; 
                ant_i < Pop[colony_i].MyAnts.size(); ++ant_i)  //actives may quit, idle may work
        {
            //if ant active 
            if (Pop[colony_i].MyAnts[ant_i].curr_act < Par.tasks)
            {
                //record last act
                Pop[colony_i].MyAnts[ant_i].last_act = 
                    Pop[colony_i].MyAnts[ant_i].curr_act; 
                
                // wanna quit?
                QuitTask(Pop[colony_i], Pop[colony_i].MyAnts[ant_i], Pop[colony_i].MyAnts[ant_i].curr_act, Par); 
            }
    else TaskChoice(Par, Pop[colony_i], Pop[colony_i].MyAnts[ant_i]); //if inactive, choose a task 

    //update number of switches after choosing tasks
    if(Pop[colony_i].MyAnts[ant_i].curr_act!=2)//if ant active 
        {
        Pop[colony_i].numacts[Pop[colony_i].MyAnts[ant_i].curr_act] += 1;
        Pop[colony_i].MyAnts[ant_i].countacts[Pop[colony_i].MyAnts[ant_i].curr_act] += 1;
        Pop[colony_i].workfor[Pop[colony_i].MyAnts[ant_i].curr_act] += Par.alfa[Pop[colony_i].MyAnts[ant_i].curr_act]; // update the work done for that task 
        if (Pop[colony_i].MyAnts[ant_i].last_act!=7 && Pop[colony_i].MyAnts[ant_i].last_act != Pop[colony_i].MyAnts[ant_i].curr_act)
            Pop[colony_i].MyAnts[ant_i].switches+=1;
        }

    /*        
    switch (Pop[colony_i].MyAnts[ant_i].curr_act)
        {
        case 0: if (Pop[colony_i].MyAnts[ant_i].last_act == 1) Pop[colony_i].MyAnts[ant_i].switches += 1; break;

        case 1: if (Pop[colony_i].MyAnts[ant_i].last_act==0) Pop[colony_i].MyAnts[ant_i].switches +=1; break; 
        
        case 2: break;
        }		
    */
    }
    }
}  // end of UpdateAnts()
//------------------------------------------------------------------------------

// update the stimuli
void UpdateStim(Population & Pop, Params & Par)   
{
    // go through all colonies
    for (unsigned int i = 0; i < Pop.size(); ++i)
    {
        // update the stimulus for each task
        for (int task = 0; task < Par.tasks; ++task)
        {
            // calculate new stimulus level
            Pop[i].newstim[task] = Pop[i].stim[task] 
                + Par.delta[task];

#ifdef SIMULTANEOUS_UPDATE
            Pop[i].newstim[task] -= (Pop[i].workfor[task]/Par.N);
#endif

            // update the stimulus
            Pop[i].stim[task] = Pop[i].newstim[task];

            if (Pop[i].stim[task] < 0)
            {
                Pop[i].stim[task] = 0;
            }
        }
    }
}

//------------------------------------------------------------------------------

void Calc_F(Population & Pop, Params & Par)
  {
  double C;
  double totacts;
  for (unsigned int i = 0; i<Pop.size(); i++)
     {
     for (unsigned int j = 0; j<Pop[i].MyAnts.size(); j++)
        {
        totacts = 0;
        C = 0;
       // cout << "Switches for ant " << j << "\t" << Pop[i].MyAnts[j].switches << endl;
      //cout << "Workperiods "<<  Pop[i].MyAnts[j].workperiods<<endl;
        //cout << "F is " << Pop[i].MyAnts[j].F << endl;
        for (int task=0; task<Par.tasks; task++)
            totacts += Pop[i].MyAnts[j].countacts[task] ;
        //cout << "number acts" << totacts << endl;
        if (totacts > 0) 
          {
          C = double(Pop[i].MyAnts[j].switches) / Pop[i].MyAnts[j].workperiods;
          Pop[i].MyAnts[j].F = 1 - 2*C;
          Pop[i].MyAnts[j].F_franjo = 1-C;
          //cout << "F is " << Pop[i].MyAnts[j].F << endl;
          }
        }
          
     double sumF=0;
     double sumF_franjo=0;
     double activ=0;
     for (unsigned int j = 0; j<Pop[i].MyAnts.size(); j++)
        {
        totacts=0;
        for (int task=0; task<Par.tasks; task++)
            totacts += Pop[i].MyAnts[j].countacts[task] ;
        if (totacts>0 && Pop[i].MyAnts[j].curr_act!=2)
                {
                sumF += Pop[i].MyAnts[j].F;
            
                activ +=1;
                sumF_franjo += Pop[i].MyAnts[j].F_franjo;
                }
        }
    // cout << "sumF= " <<sumF<< "\t"<< "active workers= " << activ << endl; 
     Pop[i].mean_F = sumF/activ;
     Pop[i].mean_F_franjo = sumF_franjo/activ; 

     }
  }
//------------------------------------------------------------------------------

void CalcFitness(Population & Pop, Params & Par)
  {
  sum_Fit = 0;
  double total = 0;

  for (unsigned int i = 0; i<Pop.size(); i++)
     {
     Pop[i].idle = 0;
     total = Pop[i].last_half_acts[0] + Pop[i].last_half_acts[1];
     if (total == 0) Pop[i].fitness =0;
     else if (total > 0)
       Pop[i].fitness =  total * (pow((Pop[i].last_half_acts[0]/total), Par.beta_fit) * pow((Pop[i].last_half_acts[1]/total), Par.gama_fit));

     for (unsigned int j =0; j< Pop[i].MyAnts.size(); j++)
       {
       if (Pop[i].MyAnts[j].countacts [0] == 0 && Pop[i].MyAnts[j].countacts [1] == 0)
         Pop[i].idle +=1;
       }
     }

     int min_fit = Pop[0].ID;
     for (unsigned int col=1; col< Pop.size(); col++)
        {
        if (Pop[col].fitness < Pop[min_fit].fitness)
                min_fit = col;
        }

     for (unsigned int i = 0; i<Pop.size(); i++)
	{
        //cout << i <<"\t"<<  Pop[i].fitness << endl;
        Pop[i].diff_fit = Pop[i].fitness - Pop[min_fit].fitness;
        sum_Fit += Pop[i].diff_fit;
        }
        //cout << "min fit " << min_fit << "\t"<< Pop[min_fit].fitness << endl;
        //getch();

     for (unsigned int i = 0; i<Pop.size(); i++)
        {
	Pop[i].rel_fit = Pop[i].diff_fit / sum_Fit;
	if (i==0)
		Pop[i].cum_fit = Pop[i].rel_fit;
	else
		Pop[i].cum_fit = Pop[i-1].cum_fit + Pop[i].rel_fit;

    //    cout << "Colony " << i << " cum_fit = " << Pop[i].cum_fit << endl;
	}
 // assert(Pop.back().cum_fit==1);

  } // end of CalcFitness()
//-------------------------------------------------------------------------------

void Categorize(Population & Pop)
  {
  for (unsigned int i=0; i<Pop.size(); i++)
    {//creating bins
    Pop[i].LowF[0] = -1;
    Pop[i].HighF[0] = -0.9;
    Pop[i].CategF[0] = 0;
    Pop[i].LowT[0] = 0;
    Pop[i].HighT[0] = 0.5;
    Pop[i].CategT1[0] = 0;
    Pop[i].CategT2[0]=0;

    for (unsigned int index = 1; index < Pop[i].HighF.size(); index++)
      {
      Pop[i].LowF[index] = Pop[i].LowF[index-1]+0.1;
      Pop[i].HighF[index] = Pop[i].HighF[index-1]+0.1;
      Pop[i].CategF[index] = 0;

      Pop[i].LowT[index]= Pop[i].LowT[index-1]+0.5;
      Pop[i].HighT[index] = Pop[i].HighT[index-1]+0.5;
      Pop[i].CategT1[index] = 0;
      Pop[i].CategT2[index] = 0;
      }

    Pop[i].LowF[10] = 0;
    Pop[i].HighF[9] = 0;

    // relative frequencies
    if (Pop[i].MyAnts.size()-Pop[i].idle > 0)
      {
      for (unsigned int worker = 0; worker < Pop[i].MyAnts.size(); worker++)
        {
       for (unsigned int index = 0; index < Pop[i].CategF.size(); index++)
         {  // specialization
         if (Pop[i].MyAnts[worker].F >= Pop[i].LowF[index]
            && Pop[i].MyAnts[worker].F < Pop[i].HighF[index])
            Pop[i].CategF[index] += 1/ ((Pop[i].MyAnts.size())-(Pop[i].idle));

         //threshold 1
     //    cout << Pop[i].MyAnts[worker].threshold[0] << endl;
         if (Pop[i].MyAnts[worker].threshold[0] >= Pop[i].LowT[index]
             && Pop[i].MyAnts[worker].threshold[0] < Pop[i].HighT[index])
            {
             Pop[i].CategT1[index] += 1;
            }
         // threshold 2
    //     cout << Pop[i].MyAnts[worker].threshold[1] << endl;
         if (Pop[i].MyAnts[worker].threshold[1] >= Pop[i].LowT[index]
             && Pop[i].MyAnts[worker].threshold[1] < Pop[i].HighT[index])
            {
             Pop[i].CategT2[index] += 1;
            }

         }   // for category
        } // for worker
      } // end if




  double sum=0;
    for (unsigned int index = 0; index < Pop[i].CategF.size(); index++)
      {
      sum +=Pop[i].CategT1[index];
     // cout << Pop[i].LowT[index] << "\t" << Pop[i].HighT[index] << endl;
      }
    //  cout << sum << endl;
    }// end of for colony
  } // end of CategorizeF()

//------------------------------------------------------------------------------

int drawParent(int nCol, Population & Pop)
{
const double draw = Uniform(); // random number
int cmin=-1, cmax=nCol-1;

while (cmax - cmin != 1)
	{
	int cmid = (cmax+cmin)/2 ;
	if (draw < Pop[cmid].cum_fit) cmax = cmid;
	else cmin = cmid;
	}
return cmax;
}// end of drawParent()
//----------------------------------------------------------------------------------------------------

void MakeSexuals(Population & Pop, Params & Par)
{
 mySexuals.resize(2 * Par.Col); // number of sexuals needed
 parentCol.resize(mySexuals.size());

 for (unsigned int ind = 0; ind < mySexuals.size(); ind ++)
	{
        //initialize sexuals
        mySexuals[ind].threshold.resize(Par.tasks);
        //mySexuals[ind].act.resize(NULL);
        mySexuals[ind].countacts.resize(NULL);
        mySexuals[ind].last_act = NULL;
        mySexuals[ind].curr_act = NULL;
        mySexuals[ind].switches = NULL;
        mySexuals[ind].F =NULL;
        mySexuals[ind].mated = false;
        // draw a parent colony for each sexual
	parentCol[ind] = drawParent(Pop.size(), Pop);
      //	cout << "Sexual " << ind << " is by col " << parentCol[ind] << endl;

	Inherit(mySexuals[ind], Pop[parentCol[ind]].queen, Pop[parentCol[ind]].male, Par);

	}
}
//-------------------------------------------------------------------------------------------
void MakeColonies(Population &Pop, Params &Par)
{
int mother, father;

for (unsigned int col = 0; col < Pop.size(); col++)
	{
	do
	{
	mother = RandomNumber(mySexuals.size());
	father = RandomNumber(mySexuals.size()) ;
	}	while
		(mother == father || mySexuals[mother].mated==true || mySexuals[father].mated==true);
	mySexuals[mother].mated = true;
	mySexuals[father].mated = true;

	Pop[col].queen = mySexuals[mother]; Pop[col].male = mySexuals[father]; // copy new males and females

	} // end for Colonies


} // end MakeColonies()
//-----------------------------------------------------------------------------------------------------

bool compare_vector_elements(vector<int>anyvector, int size_vector)
{
 bool mybool;   
 for (int m = 0; m < size_vector-1; m++) 
            {
                for (int n = m+1; n<size_vector; n++)
                    {
                    if (anyvector[m] == anyvector[n])
                        {
                        mybool=false;
                    
                   //     cout << anyvector[m] << " = " << anyvector[n] << endl;
                        break;
                        }
                    else mybool=true;
                    }
            }

return mybool;
}

//------------------------------------------------------------------------------------------------------


int main(int argc, char* argv[])
{
    Params myPars;

    ifstream inp("params.txt");

    myPars.InitParams(inp);

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
    
    out << "Gen" << "\t" << "Col"  << "\t" << "NumActs1" << "\t" << "NumActs2" <<
	"\t" << "WorkAlloc1" << "\t" << "WorkAlloc2" <<"\t" << "Idle"<< "\t" << "Fitness" <<endl; 

    out3 << "Gen" << ";" << "Time" << ";" << "Col" << ";" << "Stim1" << ";" << "Stim2" << ";" << "Workers1" << ";" << "Workers2" << ";" << "Fitness" << endl;
    }
    
    
else 
    {
    out4.open("data_1gen.txt");

    out4 << "Time" << ";" << "Col" << ";" << "Stim1" << ";" << "Stim2" << ";" << "Workers1" << ";" << "Workers2" << ";" << "Fitness" << ";" 
         << "Mean_F" << ";" << "Mean_F_franjo" << endl;
    } 

out2.open("threshold_distribution.txt");
    
    
SetSeed(myPars.seed);

//headers for data
Population MyColonies;
InitFounders(MyColonies, myPars);

for (int g= 0; g<myPars.maxgen; g++)
{
//cout << "Generation " << g << endl;

Init(MyColonies, myPars);
//cout<< "Colonies initialized" << endl;
double equil_steps=0;

// pick ten random colonies

vector<int>random_colonies(10);

#ifdef DEBUG
cout << MyColonies.size() <<endl;
cout <<myPars.maxtime << endl;

#endif
/*
vector<int>colnumbers(MyColonies.size());
for (unsigned col =0; col<MyColonies.size(); col++)
    colnumbers[col]=col;

if (g>= 40 && g<=60)
	{
    random_shuffle(colnumbers.begin(),colnumbers.end());
    int start = RandomNumber(colnumbers.size()-random_colonies.size());
    
	for (unsigned int rcol = 0; rcol < random_colonies.size(); rcol++)        
		    random_colonies[rcol] = colnumbers[start+rcol];
//	for (unsigned int rcol = 0; rcol < random_colonies.size(); rcol++)
//		    cout << "Random colonies " << random_colonies[rcol] << endl;
    
    }
*/

for (int k = 0; k < myPars.maxtime; k++)
    {
    //cout << "Time step" << k << "\t";

    UpdateAnts(MyColonies, myPars);
    /*
    if (g>=40 && g<=60 && myPars.Col>1 )
	    {
	    for (unsigned int rcol = 0; rcol < random_colonies.size(); rcol++)
		    {
		    out3 << g << ";" << k << ";" << random_colonies[rcol] << ";" ; 
		     for (int task=0; task<myPars.tasks; task++)
                out3 << MyColonies[random_colonies[rcol]].stim[task] << ";"; 
		     for (int task=0; task<myPars.tasks; task++)
	            out3 << MyColonies[random_colonies[rcol]].workfor[task]/myPars.alfa[task] << ";";
             out3 << MyColonies[random_colonies[rcol]].fitness << endl;	     
		    }
	    }
    */


    UpdateStim(MyColonies, myPars);
   
    if (k >= myPars.maxtime/2)
    	{
		for (unsigned int col = 0; col < MyColonies.size(); col++)
			{
            for (int task=0; task<myPars.tasks; task++) 
                MyColonies[col].last_half_acts[task] += MyColonies[col].workfor[task]/myPars.alfa[task];
			}	
        }
	
    
    if(k >= myPars.maxtime-51) //calculate work allocation for the last fifty timesteps
        {
		equil_steps +=1;
         for (unsigned int col = 0; col < MyColonies.size(); col++)
	        {
		    for (int task=0; task<myPars.tasks; task++)
                MyColonies[col].mean_work_alloc[task]+=MyColonies[col].workfor[task]/myPars.alfa[task];
	        }

	    }
	
    
    Calc_F(MyColonies, myPars);
    if (k == myPars.maxtime-1)
         {
         CalcFitness(MyColonies, myPars);
        // Categorize(MyColonies);
	  
         
         if ((g<=100 || g%10==0) && myPars.Col>1)
            {
            //cout <<  "Generation " << g << endl;
            for (unsigned int col = 0; col < MyColonies.size(); col++)
                {

		        for (int task=0; task<myPars.tasks; task++)
                    MyColonies[col].mean_work_alloc[task]/=equil_steps;
                    
                out << g << "\t" << col << "\t"; 
		        for (int task=0; task<myPars.tasks; task++)
                    out << MyColonies[col].numacts[task] << "\t";

		        for (int task=0; task<myPars.tasks; task++)
                    out << MyColonies[col].mean_work_alloc[task] << "\t"; 
		        
		        out << "\t" << MyColonies[col].idle << "\t" << MyColonies[col].fitness << endl; 
	           // output only thresholds of foundresses and mean specialization 
			    out2 << g << ";";
		            for (int task=0; task<myPars.tasks; task++)
				        out2 << MyColonies[col].male.threshold[task] << ";"; 
		            for (int task=0; task<myPars.tasks; task++)
				        out2 << MyColonies[col].queen.threshold[task] << ";"; 
                    out2 << MyColonies[col].mean_F <<";"<< MyColonies[col].mean_F_franjo << endl;
	    	        }
                }
            }
         

     // one colony only
   if (myPars.Col==1)
    {
    for (unsigned int col = 0; col < MyColonies.size(); col++)
            {
                out4 << k << ";" << col << ";"; 
		        for (int task=0; task<myPars.tasks; task++)
		            out4 <<MyColonies[col].stim[task] << ";"; 
		        for (int task=0; task<myPars.tasks; task++)
	                out4 << MyColonies[col].workfor[task]/myPars.alfa[task] << ";";
	            out4 << MyColonies[col].fitness << ";" << MyColonies[col].mean_F << ";" << MyColonies[col].mean_F_franjo << endl;	     
                
                for (unsigned int ind = 0; ind < MyColonies[col].MyAnts.size(); ind++)
	    	        {
			        out2 << g << ";";
		            for (int task=0; task<myPars.tasks; task++)
				        out2 << MyColonies[col].MyAnts[ind].threshold[task] << ";"; 
                    out2 << MyColonies[col].MyAnts[ind].F << ";"<< MyColonies[col].MyAnts[ind].F_franjo <<endl;
	    	        }
            }
    }

   /*  // average colonies
    double mean_work1, mean_work2;
    double mean_stim1, mean_stim2;
    double sum_work1 = 0; double sum_work2 = 0;
    double sum_stim1 = 0; double sum_stim2=0;
    for (unsigned int i = 0; i<MyColonies.size(); i++)
          {
          sum_work1 += MyColonies[i].workfor1/myPars.alfa1;
          sum_work2 += MyColonies[i].workfor2/myPars.alfa2;
          sum_stim1 += MyColonies[i].stim1;
          sum_stim2 += MyColonies[i].stim2;
          }
    mean_work1 = sum_work1/MyColonies.size();
    mean_work2 = sum_work2/MyColonies.size();
    mean_stim1 = sum_stim1/MyColonies.size();
    mean_stim2 = sum_stim2/MyColonies.size();
   */ 
   /*
    out << g << "\t" << k << "\t"  << mean_work1 << "\t"
          << mean_work2 << "\t" << mean_stim1 << "\t"
          << mean_stim2  << endl;
    */
    }

MakeSexuals(MyColonies, myPars);
MakeColonies(MyColonies, myPars);

}

//cout << "done!" << endl;

}
//---------------------------------------------------------------------------



