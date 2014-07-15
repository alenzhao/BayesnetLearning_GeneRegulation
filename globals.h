/*	globals.h

	Declarations of all common global variables.
*/


#ifndef GLOBALS_H
#define GLOBALS_H

#include "badefs.h"

using namespace std;

// **************** Some global variables *****************
extern const double func_depths[];
extern int nfunc;
extern const int tss_thrds[];
extern int ntsst;
extern const int dist_thrds[];
extern int ndistt;
extern const int loop_thrds[];
extern int nloopt;
extern int motifcand;
extern double logK;
extern int prior;
extern set<string> primo;
extern int pricnt;
extern map<string, int> mtss;
extern string rb;
extern bool itag;
// *********** Motif binding global variables **********
extern MotifMap allbind;
extern set<int> mbnd;
extern vector<MotifScore> mscor;


#define VERBOSE	// verbose mode.
//#define YY1FUNC	// use 0.01-0.04 as functional depth.


#endif

