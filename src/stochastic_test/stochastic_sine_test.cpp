#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;


//Random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *rng_global; // gnu scientific rng 

//Definining parameters
struct Params
{
	double A; //Deterministic factor
	double B; //Stochastic factor
	int maxtime; //timesteps per generation
	int genspercycle; //Generations per environmental cycle
	int randommax; //Maximum value of positive random number
	int seed;
  double delta;
  int gensdone = 0; //Generations completed
  int stepsdone = 0; //Timesteps completed in current generation
 
 //Function to read in parameters via stream
 istream & InitParams(istream & inp);
};


// A stream that reads in the parameter file
istream & Params::InitParams(istream & in)
{
	in >> A //Deterministic factor
		>> B //Stochastic factor
		>> maxtime //timesteps per generation
		>> genspercycle //Generations per environmental cycle
		>> randommax //Maximum value of positive random number
		>> seed;
   return in;
};

// Diagnostic function: prints the results of the stream to verify that it performed correctly.
void ShowParams(Params & Par)
{
	cout << "Deterministic Factor " << Par.A << endl;
	cout << "Stochastic Factor " << Par.B << endl;
	cout << "Timesteps per Generation " << Par.maxtime << endl;
	cout << "Generations per Cycle " << Par.genspercycle << endl;
	cout << "Maximum Random Number " << Par.randommax << endl;
	cout << "Seed " << Par.seed << endl;
};

// Performs stochastic sine equation (Botero et al. 2015)
void Stochsine(Params & Par)
{
	Par.delta = Par.A  * sin((2 *

		//pi
		3.14159265358979323846 * 
    
    //Calculate cumulative timesteps
    ((Par.maxtime*Par.gensdone)+Par.stepsdone)
    
    ) / Par.maxtime* Par.genspercycle)
		+ Par.B * 
    
    //Random number between 0 and randdommax
    gsl_rng_uniform_pos(rng_global) * Par.randommax;
}

int main()
{
	//Load up the parameters and print them
	Params myPars;
	ifstream inp("env_params.txt");
	myPars.InitParams(inp);
	ShowParams(myPars);

	// set up the random number generators from the gnu gsl library
	gsl_rng_env_setup();
	T = gsl_rng_default;
 
	rng_global = gsl_rng_alloc(T);
	gsl_rng_set(rng_global, myPars.seed);

	// Run the formula
	Stochsine(myPars);

	// Print delta
		cout << "Delta: " << myPars.delta << endl;

	return 0;
}