/*	globals.cpp

	Definitions of all common global variables.
*/

#include "globals.h"

// Functional depths.
// Use 0.01-0.04 for YY1 project. Default starts from 0.05-0.95 in Beer's paper.
#ifdef YY1FUNC
const double func_depths[] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 
						0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95};
#else
const double func_depths[] = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 
						0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95};
#endif
int nfunc = sizeof func_depths/sizeof func_depths[0];

// Distance to tranlation start site.
const int tss_thrds[] = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 
						220, 240, 260, 280, 300, 320, 340, 360, 380, 400,
						420, 440, 460, 480, 500, 520, 540, 560, 580, 600};
int ntsst = sizeof tss_thrds/sizeof tss_thrds[0];

// Distance between two motifs.
const int dist_thrds[] = {20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 
						220, 240, 260, 280, 300, 320, 340, 360, 380, 400,
						420, 440, 460, 480, 500, 520, 540, 560, 580, 600};
int ndistt = sizeof dist_thrds/sizeof dist_thrds[0];

// Distance for two looping motifs.
const int loop_thrds[] = {1000, 2000, 3000, 4000, 5000};
int nloopt = sizeof loop_thrds/sizeof loop_thrds[0];

int motifcand = 50;	// maximum number of motifs to be included.
double logK = 0.0;		// parameter for network structure prior.
int prior = 0;	// flag for setting prior counts for motifs that should be included into Bayesian networks.
set<string> primo;	// map to store motifs that should be added to Bayesian network apriori.
int pricnt = 20;	// prior counts for preferred motifs.
map<string, int> mtss;	// map for translational(transcriptional) start sites.
string rb = "111110";	// rule bit-string.
bool itag = false;	// Mutual information tag.

MotifMap allbind;
set<int> mbnd;
vector<MotifScore> mscor;

