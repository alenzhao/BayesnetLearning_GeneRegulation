/*	sa.cpp

	Definitions of SA related global variables.
*/

#include "sa.h"
#include "badefs.h"

// Simulated annealing.
int Repeat = 40;	// Number of repeats. Each repeat corresponds to one temperature change.
int Iteration = 20;	// Number of iterations. Each iteration corresponds to a traverse of all candidate motifs.
int Changes = 500;	// Number of requried changes. If this is met, procedure goes to the next repeat.
double Alpha = 0.9;	// Parameter to control the rate of temperature change.
double Initemp = 10.0;	// Initial temperature.
double Temp = Initemp;	// Current temperature to control the jumping rate.
int chng = 0;	// counter for network structure changes.
int maxpa = 5;	// Maximum number of parents.
int DeterNum = 5;	// Number of genes satisfying constraints to determine deteriorate condition.
double DeterScor = 0;	// Difference of scores between current and best to determine deteriorate condition.
bool Restag = false;	// Tag to determine when to restart SA if bad condition happens.
int Resthrld = 200;	// Threshold of changes to set "Restag".
int Restarts = 10;	// Maximum number of restarts for each temperature.
int rests = 0;	// counter for restarts at each temperature. 

// Best solution.
struct BSolu bsolu;	// A structure to store the best solution found.
bool tagbests = false;	// Do not use best solution as default.
