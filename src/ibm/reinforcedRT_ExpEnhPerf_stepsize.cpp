// C++ code for the evolutionary reinforced threshold model, 
//
// see Duarte (2012) Evolution of threshold reinforcement leading to division
// of labor. Chapter 4 of PhD Thesis: 
// Evolution of Self-organized Division of Labor in Social Insects 
// University of Groningen
// 
//
// based on the Fixed threshold model
// with stimulus update per ant, or per timestep (see define)
//
//
// TODO: 
// - multiple tasks
// - make sure want_task is indeed only one task (see 

//---------------------------------------------------------------------------
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstring>
#include <termios.h>
#include <unistd.h>
#include <sys/stat.h>

//#define DEBUG
//#define SIMULTANEOUS_UPDATE
//#define STOPCODE
//#define WRITE_LASTGEN_PERSTEP
//---------------------------------------------------------------------------

using namespace std;

// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *rng_r; // gnu scientific rng 


struct Params
{
    int N; //number of workers
    int Col; // number of colonies
    int maxtime; // time steps
    double p;  // quitting probability
    int tasks;
    double mutp;//mutation probability
    int maxgen;
    double beta_fit, gamma_fit; // exponents of the fitness function
    int seed;
    double recomb;
    int timecost;
    double mutstep;
    double threshold_noise; // noise when comparing thresholds to stimulus levels
    double initStim; // initial stimulus value
    double p_wait; // probability that ant has to wait c time steps before switching
    int tau; // time step from which fitness is counted
    double initForget; // initial Forget value
    double initLearn; // initial Learn
    double step_gain_exp; // stepsize when gaining experience points
    double step_lose_exp; // stepsize when losing experience points
    double K; // speed with which efficiency increases with experience

    vector<double> meanT; // mean threshold 
    vector<double> delta; // rate of increase in stimulus 
    vector<double> alfa_max; // maximum efficiency with which work is done (when fully experienced)
    vector<double>alfa_min; // minimum efficiency with which work is done (when inexperied)
    vector<double> beta;
    
    istream & InitParams(istream & inp);
};

struct Ant
{
    // genome
    double learn;
    double forget;
    // behaviour
    //vector < bool > act; // active or not at task 1 or 2
    vector < double > threshold;
    vector < double > alfa; // Strenght with which experience level affects efficiency
    vector < int > countacts;   // counter of acts done by this ant
    vector < bool > want_task; // whether an individual would accept an offered task (does not mean it will do the task)
    vector < double > experience_points; // e_ij in Duarte 2012 chapter 5
    int last_act; // keep track of the last act that an individual did
    int curr_act; // keep track of current act an individual is doing
    int switches; // number transitions to a different task
    int workperiods; // number of working periods 
    double F; // specialization value
    bool mated; // only for queens, keep track of who is already mated
    double F_franjo; // Franjo's specialization value
    int count_time; //COUNTER OF TIMESTEPS TO SWITCH TASK
    int ID_ant; // individual ID of an ant
};

// declare populations of Workers and Sexuals
typedef vector < Ant > Workers;
typedef vector < Ant > Sexuals;


// ok, define a colonoy
struct Colony
{
    Workers MyAnts; // ants in the colony
    Ant male, queen; // king & queen
    int ID; // id of the colony
    
    vector<double> stim; // the different stimuli for the various tasks in the colony 
    vector<double> newstim; // I don't know

    vector<double> workfor;  // number acts * eff each time step
    vector < int > numacts_step; // number of acts performed per task each time step
    vector < int > numacts_total; // number of total acts performed per task

    double idle; // proportion workers that _never_ worked in the simulation 
    double inactive; // proportion workers that were idle each time step

    vector <double>fitness_work; //number of acts * eff performed in the time steps counting for fitness 
    double fitness;
    double diff_fit;// fitness difference to minimal fitness
    double rel_fit; // fitness relative to whole population
    double cum_fit; //cumulative fitness
    vector<double>mean_work_alloc;
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

// define a population of colonies
typedef vector < Colony > Population;

// Sexual individuals that are going to found a new colony
Sexuals mySexuals;

// keep track of the ID of the parental colony
vector <int> parentCol;

// some stats
double sum_Fit = 0;
int simstart_generation = 0;

// if one big simulation is broken up into several 'parts' (e.g., because it takes very long)
// denote the current part
int simpart; 


int mygetch(void)
{
    struct termios oldt, newt;
    int ch;
    tcgetattr( STDIN_FILENO, &oldt );
    newt = oldt;
    newt.c_lflag &= ~( ICANON | ECHO );
    tcsetattr( STDIN_FILENO, TCSANOW, &newt );
    ch = getchar();
    tcsetattr( STDIN_FILENO, TCSANOW, &oldt );
    return ch;
}

// initialize the parameters from a textfile file in the local folder
// which is all generated through python
istream & Params::InitParams(istream & in)
{
    tasks = 2; 
    meanT.reserve(tasks);
    delta.reserve(tasks);
    alfa_max.reserve(tasks);
    alfa_min.reserve(tasks);
    beta.reserve(tasks);

    // read in parameter values from the input stream (file)
    in >> N >> // 
        Col >> // number of colonies
        maxtime;  // number of timesteps work is performed before reproduction
   
    double tmp;

    // get initial threshold values for each task
    for (int i = 0; i < tasks; ++i) 
    {
        in >> tmp; 
        meanT.push_back(tmp);
    }

    // get increase in stimulus intensity for each task
    for (int i = 0; i < tasks; ++i) 
    {    
        in >> tmp; 
        delta.push_back(tmp);
    }

    // get maximum efficiency of work for each task
    for (int i = 0; i < tasks; ++i)
    {
        in >> tmp; 
        alfa_max.push_back(tmp);
    }

    // get minimum efficiency of work for each task
    for (int i = 0; i < tasks; ++i)
    {
        in >> tmp; 
        alfa_min.push_back(tmp);
    }

    // get stimulus decay for each task
    for (int i = 0; i < tasks; ++i) 
    {
        in >> tmp;
        beta.push_back(tmp);
    }

    in >> p >> // quitting probability
        mutp >> // mutation probability
        maxgen >> // get maximum number of generations
        beta_fit >> // exponent of the first task in fitness function
        gamma_fit >> // exponent of the second task in fitness function
        recomb >> // recombination rate
        timecost>> // switching cost dependent on the time the task is performed
        mutstep >> // standard deviation of mutational distribution
        initStim >> // initial level of the stimulus
        p_wait >> // probability that ant has to wait c time steps before switching
        tau >> // timestep from which fitness is counted
        initForget >>
        initLearn >>
        step_gain_exp >>
        step_lose_exp >>
        threshold_noise >>
        K >>
        seed;
         
    return in;
}

//=================================================================================================================
// Function to check if a file exists already
bool FileExists(string strFilename) 
{
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  
  return(blnReturn);
}


//===========================================================================================
//CopyFile is a simple function that copies a file from arg1 to arg2

int CopyFile(string initialFilePath, string outputFilePath)
{		
	ifstream initialFile(initialFilePath.c_str(), ios::in|ios::binary);	
	ofstream outputFile(outputFilePath.c_str(), ios::out|ios::binary);	
	//defines the size of the buffer	
	initialFile.seekg(0, ios::end);	
	long fileSize = initialFile.tellg();	
	//Requests the buffer of the predefined size	
	//As long as both the input and output files are open...	
	if(initialFile.is_open() && outputFile.is_open())	
	    {		
	    short * buffer = new short[fileSize];		
	    //Determine the file's size		
	    //Then starts from the beginning		
	    initialFile.seekg(0, ios::beg);		
	    //Then read enough of the file to fill the buffer		
	    initialFile.read((char*)buffer, fileSize);		
	    //And then write out all that was read		
	    outputFile.write((char*)buffer, fileSize);		
	    delete[] buffer;	
	    }	

	//If there were any problems with the copying process, let the user know	
	else if(!outputFile.is_open())	
	    {		
	    cout<<"I couldn't open "<<outputFilePath<<" for copying!\n";		
	    return 0;	
	    }	
	    else if(!initialFile.is_open())	
		{		
		cout<<"I couldn't open "<<initialFilePath<<" for copying!\n";		
		return 0;	
		}			

//	initialFile.close();	
	outputFile.close();	
	return 1;
} // end of CopyFile()
//=================================================================================================================
//=================================================================================================================
// Function that initializes founders from (previous) data 

void StartFromLast (istream & in, Params & Par, Population & Pop)
{
	in >> simpart 
	   >> simstart_generation; 

	for (unsigned int i = 0; i < Pop.size(); i++)
		{
			in >> Pop[i].male.learn
                >> Pop[i].male.forget
                >> Pop[i].queen.learn
                >> Pop[i].queen.forget;
		}
}

//============================================================================================

void ShowAnts(Colony & anyCol)
{
    cout << "=======================================" << endl;
    cout << "Current values of:" << endl; 
    cout << "\t" << endl;
    for (unsigned int ant=0; ant < anyCol.MyAnts.size(); ant++)
	{
	    cout << "ant " << ant << endl;
        cout << "\t" << endl;
	    cout << "count acts " << anyCol.MyAnts[ant].countacts[0] << "\t" << anyCol.MyAnts[ant].countacts[1]  << endl;
        cout << "thresholds " << anyCol.MyAnts[ant].threshold[0] << "\t" << anyCol.MyAnts[ant].threshold[1] << endl;
        cout << "effic " << anyCol.MyAnts[ant].alfa[0] << "\t" << anyCol.MyAnts[ant].alfa[1] << endl;
	    cout << "F " << anyCol.MyAnts[ant].F << endl; // specialization value
	    cout << "switches " << anyCol.MyAnts[ant].switches << endl;
	    cout << "workperiods " << anyCol.MyAnts[ant].workperiods << endl;
        cout << "\t" << endl;
	}
}

//--------------------------------------------------------------------------------
void ShowParams(Params & Par)
{
    cout << "Workers " << Par.N << endl;
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
        cout << "max effic " << task << "\t" << Par.alfa_max[task] << endl;
    }

    for (int task=0; task<Par.tasks; task++)
    {
        cout << "min effic " << task << "\t" << Par.alfa_min[task] << endl;
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
    cout << "Exp task 2 " <<  Par.gamma_fit << endl;
    cout << "seed " << Par.seed << endl;
    cout << "recombination " << Par.recomb << endl;
    cout << "timecost " << Par.timecost << endl;
    cout << "mutstep " << Par.mutstep << endl;
    cout << "tau " << Par.tau << endl;
    cout << "initial Stimulus " << Par.initStim << endl;
    cout << "p_wait " << Par.p_wait << endl;
    cout << "initial Learn " << Par.initLearn << endl;
    cout << "initial Forget " << Par.initForget << endl;
    cout << "stepsize gain exp " << Par.step_gain_exp << endl;
    cout << "stepsize to lose exp " << Par.step_lose_exp << endl;
    cout << "K" << Par.K << endl;
}
//========================================================================================
// Function to initialize founders with their own thresholds and Learn and Forget parameters
void InitFounders(Population &Pop, Params &Par)
	{
#ifdef DEBUG  
        // some bits are taylored for two tasks so we have to throw this assert for now
        assert(Par.tasks==2);
        cout <<Par.Col << endl;
#endif
        Pop.resize(Par.Col);
        //cout << "Initializing founders!" << endl;
	for (unsigned int i = 0; i < Pop.size(); i++)
		{

            // now specify the sizes of the threshold vectors for the 
            // number of tasks which you define in the parameter file
		Pop[i].male.threshold.resize(Par.tasks);
		Pop[i].queen.threshold.resize(Par.tasks);

        // set the initial thresholds for each individual
		for(int task=0; task<Par.tasks; task++)
		    {
                // option 1: all founders have equal thresholds
		    Pop[i].male.threshold[task]= Par.meanT[task];
	            Pop[i].queen.threshold[task]= Par.meanT[task];
                // option 2: founders have thresholds drawn from normal distribution
		    //Pop[i].male.threshold[task]= Normal(Par.meanT[task],1);
			    //Pop[i].queen.threshold[task]= Normal(Par.meanT[task],1);
		    }

        // set the initial learn and forget values for the colony
        Pop[i].male.learn = Par.initLearn;
        Pop[i].male.forget = Par.initForget;
        Pop[i].queen.learn = Par.initLearn;
        Pop[i].queen.forget = Par.initForget;

		Pop[i].male.mated =true;
		Pop[i].queen.mated = true;
               // cout << Pop[i].male.threshold[0] << "\t" << Pop[i].male.threshold[1] << endl;
               // cout << Pop[i].queen.threshold[0] << "\t" << Pop[i].queen.threshold[1] << endl;
                //getch();
		}
	}
//----------------------------------------------------------------------------------------------------------------------
//=============================================================================
//end of InitFounders()

void Mutation(double & trait, double & parent, Params &Par)
    {
        if (Par.mutp > gsl_rng_uniform(rng_r))
            trait = parent + gsl_ran_gaussian(rng_r, Par.mutstep);
        else trait = parent;
    }

//==============================================================================
//end of Mutation()
void Inherit(Ant &Daughter, Ant &Mom, Ant &Dad, Params &Par)
{
    // recombination rate is sampled from a uniform distribution?
    // TODO
    double rec = gsl_rng_uniform(rng_r);

    if (Par.recomb > rec)
        {

            // which of the parents whill inherit the first trait
        int who = gsl_rng_uniform_int(rng_r, 2);

           //each task will be inherited from different parents
        if (who ==0)
                {
                    Mutation(Daughter.learn, Mom.learn, Par);
                    Mutation(Daughter.forget, Dad.forget, Par); 
                }
        else 
                {
                    Mutation(Daughter.learn, Dad.learn, Par);
                    Mutation(Daughter.forget, Mom.forget, Par);
                }
    
        } 
                //assert(Daughter.threshold[task]>0);
    else 
        {
       //both tasks will be inherited from the same parent 
        int who = gsl_rng_uniform_int(rng_r, 2);
        
                //cout << "For task " << task << endl;
                //cout << anynumber << endl;
                  //  cout << "inherited from " << who << endl;
            if (who ==0)
                {
                    Mutation(Daughter.learn, Mom.learn, Par);
                    Mutation(Daughter.forget, Mom.forget, Par);
                }

            else if (who==1)
                {
                    Mutation(Daughter.learn, Dad.learn, Par);
                    Mutation(Daughter.forget, Dad.forget, Par);
                }

//            Daughter.learn = gsl_rng_uniform(rng_r) < 0.5 ? Mutate(Dad.Learn) : Mutate(Mother.Learn)
//            Daughter.forget = gsl_rng_uniform(rng_r) < 0.5 ? Mutate(Dad.Forget) : Mutate(Mother.Forget)
        }

        if (Daughter.learn < 0) Daughter.learn = 0;
        if (Daughter.forget < 0) Daughter.forget =0;

} // end of Inherit

//=======================================================================================

// According to formalas in Duarte Thesis
void UpdateEfficiency(Ant & anyAnt, Params & Par)
{
    for (int job = 0; job < Par.tasks; job ++)
        {
        double tmp_1 = Par.K * anyAnt.experience_points[job];
        double tmp_2 = Par.alfa_min[job] * exp(tmp_1);
        anyAnt.alfa[job] = Par.alfa_max[job] * tmp_2 / (tmp_2 + (1 - Par.alfa_min[job])); 
        }

}
// end UpdateEfficiency
//=======================================================================================================================
//


// now initialize an ant
void InitAnts(Ant & myAnt, Params & Par, Colony & myCol, int numID)
    {
    myAnt.ID_ant = numID;
    myAnt.threshold.resize(Par.tasks);
    myAnt.alfa.resize(Par.tasks);
    myAnt.experience_points.resize(Par.tasks);

    for (int task=0; task<Par.tasks; task++)
        {
        myAnt.threshold[task]= Par.meanT[task];
        myAnt.experience_points[task]= 0;
   //     myAnt.alfa[task] = Par.alfa_min[task];
        }

    if (Par.maxgen > 1)
        Inherit(myAnt, myCol.queen, myCol.male, Par);
    else 
        {
            myAnt.learn= Par.initLearn;
            myAnt.forget = Par.initForget;
        }

    myAnt.countacts.resize(Par.tasks);
    myAnt.want_task.resize(Par.tasks);

    for (int task = 0; task < Par.tasks; task++)
    {
        myAnt.countacts[task]=0;
        myAnt.want_task[task]=false;
    }

    myAnt.last_act = 7; // initiate it at an impossible value for a task, because 0 is a task
    myAnt.curr_act = 2; // 0 = task 1, 1 = task 2, 2 = idle
    myAnt.switches = 0;
    myAnt.workperiods=0;
    myAnt.F = 10;
    myAnt.F_franjo = 10;
    myAnt.mated = false;
    myAnt.count_time=0;
    UpdateEfficiency(myAnt, Par); 
    }


//------------------------------------------------------------------------------------

// initialize the population
void Init(Population & Pop, Params & Par)
{
#ifdef DEBUG
    // when debugging, output some stuff
    cout << Pop.size() << endl;
    cout << Par.N << endl;
    assert(Pop.size()==Par.Col);
#endif
      
    //cout << " Initiating colonies" << endl;
    for (unsigned int colony_i = 0; colony_i < Pop.size(); ++colony_i)
    {
        // resize the colony population to fit N individuals
        Pop[colony_i].MyAnts.resize(Par.N);

        // give colony particular id (for debugging purposes)
        Pop[colony_i].ID = colony_i;

        // set colony fitness to 0
        Pop[colony_i].fitness = 0;
        Pop[colony_i].rel_fit = 0;
        Pop[colony_i].cum_fit = 0;
        Pop[colony_i].diff_fit = 0;

        // set counter of idle workers to 0
        Pop[colony_i].idle = 0;
        Pop[colony_i].inactive = 0;

        // various specialization measurements
        Pop[colony_i].mean_F = 10; 
        Pop[colony_i].var_F=0;
        Pop[colony_i].mean_F_franjo = 10;
        Pop[colony_i].var_F_franjo = 0;
        Pop[colony_i].mean_switches = 0;
        Pop[colony_i].var_switches = 0;
        Pop[colony_i].mean_workperiods=0;
        Pop[colony_i].var_workperiods=0;

        // empty pools of workers, statistics, etc
        Pop[colony_i].workfor.erase(
                Pop[colony_i].workfor.begin(),
                Pop[colony_i].workfor.end());

        Pop[colony_i].fitness_work.erase(
                Pop[colony_i].fitness_work.begin(),
                Pop[colony_i].fitness_work.end());

        Pop[colony_i].numacts_step.erase(
                Pop[colony_i].numacts_step.begin(),
                Pop[colony_i].numacts_step.end());

        Pop[colony_i].numacts_total.erase(
                Pop[colony_i].numacts_total.begin(),
                Pop[colony_i].numacts_total.end());

        Pop[colony_i].stim.erase(
                Pop[colony_i].stim.begin(),
                Pop[colony_i].stim.end());

        Pop[colony_i].newstim.erase(
                Pop[colony_i].newstim.begin(),
                Pop[colony_i].newstim.end());

        Pop[colony_i].mean_work_alloc.erase(
                Pop[colony_i].mean_work_alloc.begin(),
                Pop[colony_i].mean_work_alloc.end());

        // allocate space in the various arrays
        Pop[colony_i].workfor.reserve(Par.tasks);
        Pop[colony_i].fitness_work.reserve(Par.tasks);
        Pop[colony_i].numacts_step.reserve(Par.tasks);
        Pop[colony_i].numacts_total.reserve(Par.tasks);
        Pop[colony_i].stim.reserve(Par.tasks);
        Pop[colony_i].newstim.reserve(Par.tasks);
        Pop[colony_i].mean_work_alloc.reserve(Par.tasks);

        // put an initial 0 in the vector
        for (int task_i = 0; task_i < Par.tasks; ++task_i)
        {
            Pop[colony_i].workfor.push_back(0);
            Pop[colony_i].fitness_work.push_back(0);
            Pop[colony_i].numacts_step.push_back(0);
            Pop[colony_i].numacts_total.push_back(0);
            Pop[colony_i].stim.push_back(Par.initStim);
            Pop[colony_i].newstim.push_back(0);
            Pop[colony_i].mean_work_alloc.push_back(0);
        }

        // go through all ants in the colony and initialize the
        // individual ants
        for (unsigned int ant_i = 0; 
                ant_i < Pop[colony_i].MyAnts.size(); 
                ++ant_i)
        {
            InitAnts(Pop[colony_i].MyAnts[ant_i], Par, Pop[colony_i], ant_i);
        }
    } // end of for Pop
} // end of Init()
//==================================================================================================================
//
void UpdateStimPerAnt(Params & Par, Colony & anyCol, Ant & anyAnt, int task)
{

    //cout << "Updating the stimulus per ant" << endl;
#ifdef DEBUG
   
cout << Par.N << endl;
cout << anyCol.workfor[task] << endl;
cout <<anyCol.numacts_step[task]<< endl;
cout << anyCol.stim[task] << endl;
#endif
		anyCol.workfor[task]+=anyAnt.alfa[task];   //update workdone for fitness 
		anyCol.stim[task]-=(anyAnt.alfa[task]/Par.N); //update stimulus     
    
    if(anyCol.stim[task]<0)
        anyCol.stim[task]=0;
}

//====================================================================================================================
void UpdateThresholds_And_Experience (Ant & anyAnt, Params & Par)
    {
        switch(anyAnt.curr_act)
            {
                case 0: anyAnt.threshold[0] -= anyAnt.learn ; 
                        anyAnt.threshold[1] += anyAnt.forget;
                        anyAnt.experience_points[0]+=Par.step_gain_exp;
                        anyAnt.experience_points[1]-=Par.step_lose_exp;
                       break;

                case 1: anyAnt.threshold[1] -= anyAnt.learn; 
                        anyAnt.threshold[0] += anyAnt.forget;
                        anyAnt.experience_points[0]-= Par.step_lose_exp;
                        anyAnt.experience_points[1]+=Par.step_gain_exp;
                        break;

                case 2: anyAnt.threshold[0] += anyAnt.forget;
                        anyAnt.threshold[1] += anyAnt.forget;
                         anyAnt.experience_points[0]-=Par.step_lose_exp;
                         anyAnt.experience_points[1]-=Par.step_lose_exp;
                        break;

                default: cout << "ERROR -- ANT NOT ASSIGNED TASK NOR IDLE" << endl;
            }
        for (int task = 0; task < Par.tasks; task ++)
            {
            if(anyAnt.threshold[task]<0) anyAnt.threshold[task] =0;
            if(anyAnt.experience_points[task] < 0) anyAnt.experience_points[task] = 0;
            }
    }   
//========================================================================================================================
void UpdateSwitches(Ant & anyAnt, Params & Par)
{

    if(anyAnt.curr_act != 2) // if ant active
        {
            if(anyAnt.last_act != 7 && anyAnt.last_act != anyAnt.curr_act)
                anyAnt.switches +=1;

        } 
}
//end UpdateSwitches
//=========================================================================================================================
void QuitTask(Colony & anyCol, Ant & anyAnt, int job, Params & Par)
    {
#ifdef DEBUG
        cout << "Quitting tasks" << endl;
        cout << "chance to quit: " << Par.p << endl;

#endif
        double q = gsl_rng_uniform(rng_r);

#ifdef SIMULTANEOUS_UPDATE

        if (q <= Par.p)
         {
         //anyAnt.act[job] = false;
         //anyAnt.act[2] = true; // she quits task and becomes idle
         anyAnt.want_task[anyAnt.curr_act]=false;
         anyAnt.count_time=0; // reset time to zero, she may choose the same or another task next
         anyAnt.curr_act = 2;
         //cout << "ant quits" << endl;
         }
#endif  

#ifndef SIMULTANEOUS_UPDATE  
        if (q <= Par.p)
         {
         //anyAnt.act[job] = false;
         //anyAnt.act[2] = true; // she quits task and becomes idle
         anyAnt.want_task[anyAnt.curr_act]=false;
         anyAnt.count_time=0; // reset time to zero, she may choose the same or another task next
         anyAnt.curr_act = 2;
         //cout << "ant quits" << endl;
         }
        else UpdateStimPerAnt(Par, anyCol, anyAnt, job);
#endif

    }

//------------------------------------------------------------------------------
void DoTask ( Params Par, Colony & anyCol, Ant & anyAnt,int job)
{

         anyAnt.curr_act = job; 
         anyAnt.workperiods +=1; 
         anyAnt.countacts[job] +=1;
        
#ifndef SIMULTANEOUS_UPDATE 
         UpdateStimPerAnt(Par, anyCol, anyAnt, job);
#endif
}
//end DoTask
//---------------------------------------------------------------------------------

// let ant evaluate threshold and see if she wants to perform a task
// several outcomes: ant may prefer one or multiple tasks. In the latter
// case, one of those tasks is selected as the preferred task
// she may also want to prefer no task yet
void WantTask (Params Par, 
        Colony & anyCol, 
        Ant & focalAnt)
{
    // make a list of all the task that this ants wants to do
    // and reserve space for it
    vector<int>wanted_task_ids;
    wanted_task_ids.reserve(Par.tasks); 

    // variable to store the focal ant's threshold value + noise for a task
    double t_noise;

    // loop through all tasks and calculate thresholds
    for (int job = 0; job < Par.tasks; ++job)
    {
        // calculate threshold + random noise
        t_noise = focalAnt.threshold[job] + gsl_ran_gaussian(rng_r, Par.threshold_noise);

        // threshold cannot be negative
        if (t_noise < 0)
        {
            t_noise = 0;
        }

        // ants want to work on tasks for which 
        // - their threshold exceeds the threshold + noise
        // - the stimulus level is nonzero (i.e., work needs to be done)
        if (anyCol.stim[job] >= t_noise && anyCol.stim[job] > 0) 
        {
            // store the wanted task
            wanted_task_ids.push_back(job);  
        } 
        else // ok ant does not want this task
        {
            // if ant's threshold not high enough, then quit with wanting task
            focalAnt.want_task[job] = false;
        }
    }

    // if more than one task is above threshold, select a random task
    // that ant wants to perform
    if (wanted_task_ids.size() > 1)
    {
        int job = gsl_rng_uniform_int(rng_r, wanted_task_ids.size());
        focalAnt.want_task[wanted_task_ids[job]] = true;
    }
    else if (!wanted_task_ids.empty())
    {
        assert(wanted_task_ids[0] >= 0);
        assert(wanted_task_ids[0] < Par.tasks);

        focalAnt.want_task[wanted_task_ids[0]] = true;
    }
}
//end WantTask
//-----------------------------------------------------------------------------

// evaluate whether ant should switch tasks
void EvalTaskSwitch(Params & Par, Colony & anyCol, Ant & anyAnt,int myjob)
{
    // if it was doing this job previously 
    // or it did not do anything before
    // just perform the task
    if (myjob == anyAnt.last_act || anyAnt.last_act == 7) 
    {
        DoTask(Par, anyCol, anyAnt, myjob);
    }
    else  // ok, ant <is doing a different task than dei
    {
        // with a certain probability 
        if (Par.p_wait >= gsl_rng_uniform(rng_r) 
                && anyAnt.count_time < Par.timecost)    
        {
            //cout << "Ant wants to change task!" << endl;    
            anyAnt.curr_act=2; // stays idle for as long as count_time<timecost    
            ++anyAnt.count_time;
        }
        else
        {
            DoTask(Par, anyCol, anyAnt, myjob);
        }
    }
}
//end EvalTaskSwitch
//-------------------------------------------------------------------------------

// act of choosing a task
void TaskChoice(
        Params & Par, // parameter object
        Colony & anyCol, // current colony
        Ant & focalAnt) // the ant in question
{ 
#ifdef DEBUG

    // little piece of code to assert that ants do not want to do
    // multiple tasks at the same time
    bool wants_task = false;
    for (int task_i = 0; task_i < Par.tasks; ++task_i)
    {
        if (focalAnt.want_task[task_i])
        {
            if (!wants_task)
            {
                wants_task = true;
            }
            else
            {
                cout << "error: ant wants multiple tasks simultaneously";
                exit(1);
            }
        }
    }
#endif

     // find out if ant wants to perform a task
    for (int task_i = 0; task_i < Par.tasks; ++task_i)
    {
        // yes, ant wants to perform task so let's do it
        if (focalAnt.want_task[task_i])
        {
            EvalTaskSwitch(Par, anyCol, focalAnt, task_i); 
            return;
        }
    }

    // ant does not want to perform a task

    // make ant want task
    WantTask(Par, anyCol, focalAnt);

     // find out if ant now wants to perform a task
    for (int task_i = 0; task_i < Par.tasks; ++task_i)
    {
        // yes, ant wants to perform task so let's do it
        if (focalAnt.want_task[task_i])
        {
            EvalTaskSwitch(Par, anyCol, focalAnt, task_i); 
            return;
        }
    }
} // end of TaskChoice()
//===============================================================================================

// update all the stimulus values for each colony
void UpdateAnts(Population & Pop, Params & Par)
{
    // go through all the individual colonies
    for (unsigned int i = 0; i < Pop.size(); ++i)
    {
            for (int task =0; task<Par.tasks; task++)
             {
             Pop[i].workfor[task]=0; // reset work for tasks
             Pop[i].numacts_step[task] = 0;
             }
         random_shuffle(Pop[i].MyAnts.begin(), Pop[i].MyAnts.end());
         
         for (unsigned int j = 0; j<Pop[i].MyAnts.size(); j++)  //actives may quit, idle may work
            {
            //cout << j << " of "<< Pop[i].MyAnts.size() << " workers " << endl;
            if(Pop[i].MyAnts[j].curr_act!=2)//if ant active 
                {
                //cout << "Ant busy " << endl;
                Pop[i].MyAnts[j].last_act = Pop[i].MyAnts[j].curr_act; //record last act
                QuitTask(Pop[i], Pop[i].MyAnts[j], Pop[i].MyAnts[j].curr_act, Par); // wanna quit?
                }
            if(Pop[i].MyAnts[j].curr_act==2) TaskChoice(Par, Pop[i], Pop[i].MyAnts[j]); //if inactive, choose a task 
            
            //update number of switches after choosing tasks
            UpdateSwitches(Pop[i].MyAnts[j], Par);
            UpdateThresholds_And_Experience(Pop[i].MyAnts[j], Par);
            UpdateEfficiency(Pop[i].MyAnts[j], Par);
            } 
          
        //ShowAnts(Pop[i]);
        //mygetch();
         
        }
}  // end of UpdateAnts()
//------------------------------------------------------------------------------

void Update_Col_Data(
        int step,  // current timestep
        Population & Pop, // the metapopulation
        Params & Par, // the parameters
        double & eq_steps // 
        )
{
    // go through all colonies and update act counters
    for (unsigned int col = 0; col < Pop.size(); ++col) 
    {
        Pop[col].inactive = 0;

        for (unsigned int ind = 0; ind < Pop[col].MyAnts.size(); ind ++)
        {
            // check whether ant is active
            if (Pop[col].MyAnts[ind].curr_act < Par.tasks)
            {
                // if active update act count
                Pop[col].numacts_step[Pop[col].MyAnts[ind].curr_act] += 1; 
            }
            else // ant inactive, count it
            {
                ++Pop[col].inactive;
            }
            
            Pop[col].inactive /= Pop[col].MyAnts.size(); // proportion inactive workers 
        }
        
        Pop[col].numacts_total[0] += Pop[col].numacts_step[0]; 
        Pop[col].numacts_total[1] += Pop[col].numacts_step[1]; 
    }

    // calculate fitness if within tau timesteps from the end
	if (step >= Par.tau) 
    {
        // update the number of timesteps that fitness is counted
        ++eq_steps;

        for (unsigned int col = 0; col < Pop.size(); ++col)
        {
            for (int task=0; task<Par.tasks; task++) 
            {
                // add the number of workers to the fitness tally
                Pop[col].fitness_work[task] += Pop[col].workfor[task];
                Pop[col].mean_work_alloc[task]+=Pop[col].numacts_step[task];
            }
        }	
    }
		
} // end of UpdateColony_data
//===================================================================================================
void UpdateStim(Population & Pop, Params & Par)   // stimulus changes const increase
{
    //cout << "updating stimulus" << endl;
    for (unsigned int i = 0; i<Pop.size(); i++)
    {
#ifdef SIMULTANEOUS_UPDATE
        for (int task=0; task<Par.tasks; task++)
        {
            // update the stimulus for this task
            Pop[i].newstim[task] = Pop[i].stim[task] + Par.delta[task] - 
                (Par.beta[task]*Pop[i].stim[task]) - (Pop[i].workfor[task]/Par.N);

            Pop[i].stim[task] = Pop[i].newstim[task];

            if (Pop[i].stim[task] < 0)
            {
                Pop[i].stim[task] =0;
            }
        }
#endif

#ifndef SIMULTANEOUS_UPDATE
        for (int task=0; task<Par.tasks; task++)
        {
            Pop[i].newstim[task] = Pop[i].stim[task] + 
                Par.delta[task] - (Par.beta[task]*Pop[i].stim[task]); 

            Pop[i].stim[task] = Pop[i].newstim[task];

            if (Pop[i].stim[task] < 0)
            {
                Pop[i].stim[task] =0;
            }
        }
#endif

    }
} // end UpdateStim()
//==============================================================================================

// calculate specialization value
void Calc_F(Population & Pop, Params & Par)
{
    double C;
   
    for (unsigned int i = 0; i<Pop.size(); i++)
        {
        Pop[i].mean_F=-10;
        Pop[i].mean_F_franjo=-10; 
        Pop[i].mean_switches=0; 
        Pop[i].mean_workperiods=0; 

        double p1 = (double)Pop[i].numacts_total[0] / (Pop[i].numacts_total[0] + Pop[i].numacts_total[1]);
        
        double p2 = (double)Pop[i].numacts_total[1] / (Pop[i].numacts_total[0] + Pop[i].numacts_total[1]);
        
        double denomin = p1*p1 + p2*p2;

        double sumF=0; double sumsquares_F=0;

        double sumF_franjo=0; double sumsquares_F_franjo=0;

        double sumsquares_switches=0; double sumsquares_workperiods=0;

        double activ=0;

        for (unsigned int j = 0; j<Pop[i].MyAnts.size(); j++)
            {
            assert(Pop[i].MyAnts[j].workperiods <= Par.maxtime);
            C = 0;
            if (Pop[i].MyAnts[j].workperiods > 1)
              {
              Pop[i].mean_switches += Pop[i].MyAnts[j].switches;
              Pop[i].mean_workperiods += Pop[i].MyAnts[j].workperiods;
              
              C = double(Pop[i].MyAnts[j].switches) / (Pop[i].MyAnts[j].workperiods - 1);
              
              Pop[i].MyAnts[j].F = 1 - 2*C;

// correction for proportion of tasks
              Pop[i].MyAnts[j].F_franjo = (1-C) / denomin;
              
              sumF += Pop[i].MyAnts[j].F;
              activ +=1;
              sumF_franjo +=Pop[i].MyAnts[j].F_franjo;
              }
            }


        Pop[i].mean_switches /= activ;

        Pop[i].mean_workperiods = Pop[i].mean_workperiods/Pop[i].MyAnts.size();
        //cout << "mean of switches " << Pop[i].mean_switches << endl;
        Pop[i].mean_F = sumF/activ;
        Pop[i].mean_F_franjo = sumF_franjo/activ;

        for (unsigned int j = 0; j<Pop[i].MyAnts.size(); j++)
            {
                sumsquares_workperiods += (Pop[i].MyAnts[j].workperiods - Pop[i].mean_workperiods)*(Pop[i].MyAnts[j].workperiods - Pop[i].mean_workperiods);

                if (Pop[i].MyAnts[j].workperiods > 1)
                    {
                    sumsquares_switches += ( Pop[i].MyAnts[j].switches - Pop[i].mean_switches)*( Pop[i].MyAnts[j].switches - Pop[i].mean_switches);

                    sumsquares_F += (Pop[i].MyAnts[j].F - Pop[i].mean_F)* (Pop[i].MyAnts[j].F - Pop[i].mean_F);

                    sumsquares_F_franjo += (Pop[i].MyAnts[j].F_franjo - Pop[i].mean_F_franjo)* (Pop[i].MyAnts[j].F_franjo - Pop[i].mean_F_franjo);
                    }
            }

        Pop[i].var_switches = sumsquares_switches / activ;
        Pop[i].var_workperiods = sumsquares_workperiods / Pop[i].MyAnts.size();
        Pop[i].var_F = sumsquares_F / activ;
        Pop[i].var_F_franjo = sumsquares_F_franjo / activ;

        //cout << "mean of switches " << Pop[i].mean_switches << endl;
        } // end for Colonies

  } // end of Calc_F()
//=======================================================================================================================

// determine fitness
void CalcFitness(Population & Pop, Params & Par)
{
    sum_Fit = 0;
    double total_work = 0;

    // calculate fitness for all colonies
    for (unsigned int i = 0; i < Pop.size(); ++i)
    {
        Pop[i].idle = 0;

        total_work = Pop[i].fitness_work[0] + Pop[i].fitness_work[1];
             
        if (total_work == 0)
        {
            Pop[i].fitness = 0;
        }
        else
        {
            // multiplicative fitness
            Pop[i].fitness = Pop[i].fitness_work[0] * Pop[i].fitness_work[1];
        }

        // calculate the number of idle ants 
        for (unsigned int j = 0; j < Pop[i].MyAnts.size(); ++j)
        {
            if (Pop[i].MyAnts[j].workperiods == 0)
            {
                Pop[i].idle +=1;
            }
        }
    }

    // determine the least fit colony
    int min_fit = Pop[0].ID;

    for (unsigned int col = 1; col < Pop.size(); ++col)
    {
        if (Pop[col].fitness < Pop[min_fit].fitness)
        {
            min_fit = col;
        }
    }

    // calculate difference of each colony relative to minimum fitness
    for (unsigned int i = 0; i < Pop.size(); i++)
    {
        Pop[i].diff_fit = Pop[i].fitness - Pop[min_fit].fitness;

        sum_Fit += Pop[i].diff_fit;
    }

    Pop[0].rel_fit  = Pop[0].diff_fit / sum_Fit;
    Pop[0].cum_fit = Pop[0].rel_fit;

    // normalize difference wrt the fitness sum
    for (unsigned int i = 1; i<Pop.size(); i++)
    {
        Pop[i].rel_fit = Pop[i].diff_fit / sum_Fit;

        Pop[i].cum_fit = Pop[i-1].cum_fit + Pop[i].rel_fit;
    }

} // end of CalcFitness()
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

int drawParent(int nCol, Population & Pop)
    {
    const double draw = gsl_rng_uniform(rng_r); 
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
            mySexuals[ind].mated = false;
            // draw a parent colony for each sexual
        parentCol[ind] = drawParent(Pop.size(), Pop);
          //	cout << "Sexual " << ind << " is by col " << parentCol[ind] << endl;

        Inherit(mySexuals[ind], Pop[parentCol[ind]].queen, Pop[parentCol[ind]].male, Par);

        }
} // end of MakeSexuals
//-------------------------------------------------------------------------------------------
void MakeColonies(Population &Pop)
{
    int mother, father;

    for (unsigned int col = 0; col < Pop.size(); col++)
        {
        do
            {
            mother = gsl_rng_uniform_int(rng_r, mySexuals.size());
            father = gsl_rng_uniform_int(rng_r, mySexuals.size()) ;
            }	while
                    (mother == father || mySexuals[mother].mated==true || mySexuals[father].mated==true);

        mySexuals[mother].mated = true;
        mySexuals[father].mated = true;

        Pop[col].queen = mySexuals[mother]; Pop[col].male = mySexuals[father]; // copy new males and females

        } // end for Colonies

} // end MakeColonies()
//-----------------------------------------------------------------------------------------------------
void WriteLastGen(int generation, Params & Par, Population & Pop) 
{
    //cout << "writing lastgen" << endl;
    // if last generation, write out the founders to a file 
    static ofstream lastgen;
    //static ofstream lastgen_backup;

    if (generation%10 == 0 || generation - simstart_generation == Par.maxgen -1) 
        {
//	    if(!FileExists("lastgen.txt"))
//		{
//			lastgen_backup.open("lastgen_backup.txt");
//		}

            lastgen.open("lastgen.txt");
            lastgen << simpart + 1 << endl;
            lastgen << generation + 1 << endl;
            for (unsigned int col = 0; col < Pop.size(); col++)
                {
                        lastgen << Pop[col].male.learn<< endl;	
                        lastgen << Pop[col].male.forget<< endl;	
                   
                        lastgen << Pop[col].queen.learn << endl;	
                        lastgen << Pop[col].queen.forget << endl;	
                    
                }

	    CopyFile("lastgen.txt", "lastgen_backup.txt");
	    lastgen.close();

        }
} // end of WriteLastGen
//========================================================================================================
//
// give names to each of the datafiles
void NameDataFiles(
        string &data1, 
        string &data2, 
        string &data3, 
        string &data4, 
        string &data5, 
        string &data6, 
        string &dataant)
{
    stringstream tmp;
    tmp << "data_work_alloc_" << simpart << ".txt";
    data1 = tmp.str();

    stringstream tmp2;
    tmp2 << "allele_distrib_" << simpart << ".txt";
    data2 = tmp2.str();

    stringstream tmp3;
    tmp3 << "f_dist_" << simpart << ".txt";
    data3 = tmp3.str();

    stringstream tmp4;
    tmp4 << "branch.txt";
    data4 = tmp4.str();

    stringstream tmp5;
    tmp5 << "data_1gen_" << simpart << ".txt";
    data5 = tmp5.str();

    stringstream tmp6;
    tmp6 << "thresholds_" << simpart << ".txt";
    data6 = tmp6.str();

    stringstream tmp7;
    tmp7 << "ant_beh_" << simpart << ".txt";
    dataant = tmp7.str();
}
//=====================================================================================================

// is this a continuation of a previous run, yes or no?
// if yes, read in the last generation of the previous run and start from there
// if no, just initialize everything
void Continue_Previous_Run_Yes_No(Params & Par, Population & Pop) 
{
    // if a lastgen.txt file is present in the current directory
    // this means it is a continuation of an older run 
	if (FileExists("lastgen.txt"))
    {
        ifstream inp("lastgen.txt");
        StartFromLast(inp, Par, Pop);
    }
	else 
    {
        simpart = 1; 
        simstart_generation = 0; 
    }
}
//=====================================================================================================

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

//-------------------------------------------------------------------------------------------------------
void Write_Col_Data(
        ofstream & mydata,
        Params & Par,
        Population & Pop,
        int gen,
        int colony)
{
    //cout << "write col data " << endl;
	mydata << gen << "\t" << colony << "\t";
	for (int task=0; task<Par.tasks; task++)
		    mydata << Pop[colony].fitness_work[task] << "\t";
	for (int task=0; task<Par.tasks; task++)
		    mydata << Pop[colony].mean_work_alloc[task] << "\t"; 				
	        mydata << Pop[colony].idle 
                << "\t" << Pop[colony].inactive 
                << "\t" << Pop[colony].fitness 
                << "\t" << Pop[colony].stim[0] 
                << "\t" << Pop[colony].stim[1] 
                << "\t" << Pop[colony].mean_switches 
                << "\t" << Pop[colony].mean_workperiods   
	    << endl;  
}
//------------------------------------------------------------------------------------------------------
// write out all the alleles to get an overview of
// the amount of within and between colony genetic variation
void Write_Alleles_Spec(
        ofstream & data_reinforcement,
        ofstream & data_f,
        Params & Par,
        Population & Pop,
        int gen,
        int colony)
{
    //cout << "write alleles and spec " << endl;
	 // output only thresholds of foundresses and mean specialization 
    data_reinforcement << gen << ";";
    data_reinforcement << Pop[colony].male.learn << ";"; 
    data_reinforcement << Pop[colony].male.forget << endl; 

	data_reinforcement << gen <<";";
    data_reinforcement << Pop[colony].queen.learn << ";"; 
    data_reinforcement << Pop[colony].queen.forget << endl; 

	data_f << gen <<";" 
                << Pop[colony].mean_F_franjo << ";" 
                << Pop[colony].mean_switches << ";" 
                << Pop[colony].mean_workperiods << ";"
                << Pop[colony].var_F_franjo << ";"
                << Pop[colony].var_switches <<";" 
                << Pop[colony].var_workperiods << endl;
} 

//====================================================================================================
void Write_Branching(ofstream & afile,
        Params & Par,
        int yesno)
{
	afile << Par.mutp << ";" 
		<< Par.mutstep << ";" 
		<< Par.timecost << ";" 
		<< yesno << endl;	
}	
//====================================================================================================
void WriteBranchFile(Population & Pop,
        int generation, 
        Params & Par,
        ofstream &file4,
        string filename)
	{ 
	file4.open(filename.c_str());
	if (Check_Spec(Pop))
		{
		int branch = 1;
		Write_Branching(file4, Par, branch);	
		}

	 else
		{
		
		int branch = 0;	
		Write_Branching(file4, Par, branch);	
		}
	file4.close();

	} // end WriteBranchFile


//==============================================================================================================================================

// add headers to the data files
void Header_data(ofstream & header, ofstream & header2)
{
	    header << "Gen" << "\t" 
		<< "Col"  << "\t" 
		<< "FitWork1" << "\t" 
		<< "FitWork2" << "\t" 
		<< "WorkAlloc1" << "\t" 
		<< "WorkAlloc2" <<"\t" 
		<< "Idle"<< "\t" 
		<< "Inactive" << "\t"
		<<"Fitness" << "\t" 
		<< "End_stim1" << "\t" 
		<< "End_stim2" << "\t"
		<< "mean_switches" << "\t"
		<< "mean_workperiods" << endl; 

#ifdef WRITE_LASTGEN_PERSTEP
        header2 << "Time" << ";" 
		<< "Col" << ";" 
		<< "Stim1" << ";" 
		<< "Stim2" << ";" 
		<< "Workers1" << ";"
		<< "Workers2" << ";" 
		<< "Fitness" << ";" 
		<< "Mean_F_franjo" << endl;
#endif   
}
//=======================================================================================================================
//
//
//writing data of last generation step by step!
void WriteData_1Gen(ofstream & mydata, 
        Population & Pop, 
        Params & Par, 
        int timestep)
{
    // loop through all colonies and write stats
    for (unsigned int col = 0; col < Pop.size(); ++col)
    {
        // print timestep and colony number
        mydata << timestep << ";" << col << ";"; 

        // plot the perceived stimulus levels per task
        for (unsigned int task=0; task < Par.tasks; ++task) 
        {
            mydata << Pop[col].stim[task] << ";"; 
        }

        // plot the number of acts for each task
        for (unsigned int task=0; task<Par.tasks; task++) 
        {
            mydata << Pop[col].numacts_step[task] << ";";
        }

        // write down colony fitness and specialization values
        mydata << Pop[col].fitness << ";" 
            << Pop[col].mean_F_franjo << endl; 
    }
}
//==========================================================================================================================

// write down individual ants
void WriteAntsBeh(Population & Pop, ofstream & mydata) 
{
    for (unsigned int col = 0; col < Pop.size(); ++col)
    {
        for (unsigned int ant = 0; ant < Pop[col].MyAnts.size(); ++ant)
        {
            mydata << col << ";" << ant << ";" 
                << Pop[col].MyAnts[ant].threshold[0] << ";" 
                << Pop[col].MyAnts[ant].threshold[1] << ";" 
                << Pop[col].MyAnts[ant].countacts[0] << ";" 
                << Pop[col].MyAnts[ant].countacts[1] << ";"
                << Pop[col].MyAnts[ant].experience_points[0] << ";"
                << Pop[col].MyAnts[ant].experience_points[1] << ";"
                << Pop[col].MyAnts[ant].alfa[0] << ";"
                << Pop[col].MyAnts[ant].alfa[1] << ";" 
                << Pop[col].MyAnts[ant].switches << ";"
                << Pop[col].MyAnts[ant].workperiods << endl;  
        }
    }
}
//=========================================================================================================
// writing ants' thresholds 
void WriteAntsThresholds(Population & Pop, ofstream & mydata, int timestep, int gen) 
{
    for (unsigned int col = 0; col < Pop.size(); ++col)
    {
        for (unsigned int ant = 0; ant < Pop[col].MyAnts.size(); ++ant)
        {
            mydata << gen << ";" << timestep << ";" << col << ";" << Pop[col].MyAnts[ant].ID_ant << ";" 
                << Pop[col].MyAnts[ant].threshold[0] << ";" 
                << Pop[col].MyAnts[ant].threshold[1] << ";" 
                << Pop[col].MyAnts[ant].countacts[0] << ";"
                << Pop[col].MyAnts[ant].countacts[1] << ";"
                << Pop[col].MyAnts[ant].experience_points[0] << ";"
                << Pop[col].MyAnts[ant].experience_points[1] << ";"
                << Pop[col].MyAnts[ant].alfa[0] << ";"
                << Pop[col].MyAnts[ant].alfa[1] << endl; 
        }
    }
}
//================================================================================
int main(int argc, char* argv[])
{
        // initialize object to store all parameters
        Params myPars;
        
        // get parameters from file
        ifstream inp("params.txt");

        // add these parameters to parameter object
        myPars.InitParams(inp);

        // set up the random number generators
        // (from the gnu gsl library)
        gsl_rng_env_setup();
        T = gsl_rng_default;
        rng_r = gsl_rng_alloc(T);
        gsl_rng_set(rng_r, myPars.seed);

        // initialize the founders of all the colonies
        Population MyColonies;
        InitFounders(MyColonies, myPars);

        // this simulation run might a a continuation of a previous 
        // simulation, for example when that simulation was broken off
        // prematurely. This function checks whether lastgen.txt (the output 
        // of the previous simulation is present and initializes the simulation
        // accordingly
        Continue_Previous_Run_Yes_No(myPars, MyColonies);

        // all the files to which data is written to
        string datafile1, 
               datafile2, 
               datafile3, 
               datafile4, 
               datafile5, 
               datafile6, 
               dataants;

        // function to give the datafiles particular names
        NameDataFiles(
                datafile1, 
                datafile2, 
                datafile3, 
                datafile4, 
                datafile5, 
                datafile6, 
                dataants);

        // the corresponding output files
        static ofstream out1; 
        static ofstream out2;
        static ofstream out3;
        static ofstream out4;
        static ofstream out5;
        static ofstream out6;
        static ofstream header1;
        static ofstream header2;
        static ofstream out_ants;

        // more than one colony, so..
        if (myPars.Col > 1) 
        {
            out1.open(datafile1.c_str());

            header1.open("header_1.txt");

#ifdef WRITE_LASTGEN_PERSTEP 
            header2.open("header2.txt");
#endif

            Header_data(header1, header2);
        }
            
        out2.open(datafile2.c_str());
        out3.open(datafile3.c_str());    

#ifdef WRITE_LASTGEN_PERSTEP 
        out_ants.open(dataants.c_str()); 
        out5.open(datafile5.c_str());
        out6.open(datafile6.c_str());
#endif

        // calculate maximum number of generations
        int maxgen = simstart_generation + myPars.maxgen;

        // now go evolve
        for (int g = simstart_generation; g < maxgen; ++g)
        {
            cout << g << endl;
            Init(MyColonies, myPars);
           
            // number of timesteps that fitness is counted
            // (e.g., when time > tau)
            double equil_steps=0;

#ifdef DEBUG
            cout << MyColonies.size() <<endl;
            cout <<myPars.maxtime << endl;

#endif

            // timesteps during colony development
            for (int k = 0; k < myPars.maxtime; ++k)
            {
                // update all the stimuli of the ants and what they are doing
                UpdateAnts(MyColonies, myPars);

                // calculate specialization values
                Calc_F(MyColonies, myPars); 

                // update statistics and if beyond tau, fitness values
                Update_Col_Data(k, MyColonies, myPars, equil_steps);	

                // calculate at the end of the timestep: the ants have done something
                // which has consequences for stimulus levels, which you update here
                UpdateStim(MyColonies, myPars);
          
                // in the last timestep
                if (k == myPars.maxtime-1)
                {
                    // calculate fitness values
                    CalcFitness(MyColonies, myPars);

                    // write everything down every 100th timestep
                    // calculating fitness is only relevant when there are
                    // multiple colonies
                    if ((g <= 100 || g % 100 == 0) && myPars.Col>1)
                    {
                        //cout <<  "Generation " << g << endl;
                        for (unsigned int col = 0; col < MyColonies.size(); ++col)
                        {
                            // take averages over all tasks
                            for (int task = 0; task < myPars.tasks; ++task)
                            {
                                MyColonies[col].mean_work_alloc[task]/=equil_steps;
                            }

                            // write out the data 
                            Write_Col_Data(out1, myPars, MyColonies, g, col);

                            // write out alleles
                            Write_Alleles_Spec(out2, 
                                    out3, 
                                    myPars, 
                                    MyColonies, 
                                    g, 
                                    col);
                        } // end of for colonies
                     } // end if generations are right
                } // end if k=maxtime

            //do you want to write out the last generation step by step?
#ifdef WRITE_LASTGEN_PERSTEP
            
                if (g == simstart_generation+myPars.maxgen-1) 
                {
                    //cout << " writing last generation data " << endl;
                    WriteData_1Gen(out5, MyColonies, myPars, k);
        //             WriteAntsThresholds(MyColonies, out6, k, g);
                    if(k == myPars.maxtime -1)
                        {	
                        WriteAntsBeh(MyColonies, out_ants);	
                        }
                    //cout << "done with writing last generation data " << endl;
                }
#endif
            } // end of for k time steps
    
    // one colony only
    WriteLastGen(g, myPars, MyColonies);

    WriteBranchFile(MyColonies, g, myPars, out4, datafile4);
        if(g < myPars.maxgen -1)
            {
            MakeSexuals(MyColonies, myPars);

            MakeColonies(MyColonies);
            }
        } // end for generations
}
