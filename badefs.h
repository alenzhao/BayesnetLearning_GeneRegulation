/*	badefs.h

	Definitions of types used in Bayesian networks.
*/


#ifndef BADEFS_H
#define BADEFS_H

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>

using namespace std;

// Struct to store each single motif's name, presence node maximum score
// and corresponding functional depth.
typedef struct{
	string name;	// motif's name.
	double score;	// maximum score when as a single presence node.
	double depth;	// current chosen functional depth.
} MotifScore;

// Information of a single binding site.
typedef struct{
	char orien;		// Orientation.
	double score;	// Matrix score.
	int loc;		// Binding location.
} GBinding;

// A vector to contain all binding sites of a gene.
struct VGB{
	vector<GBinding> e;
};

// A map to contain all genes' binding.
struct GBMap{
	map<string, VGB> e;
};

// Description of a constraint.
typedef struct{
	string desc;	// a string ID for the constraint.
	int motif0;		// index to motif0 for the constraint.
	int motif1;		// index to motif1 for the constraint.
	int para;		// parameter of the constraint.
} Constraint;

// Conditional probability table: each cell represents a parent's state and a child's state.
typedef struct{
	int k0;	// number of cases when child = 0;
	int k1;	// number of cases when child = 1;
} CPTRow;

// Case in database: gene name and its label.
typedef struct{
	string name;
	int label;
} Case;

// A map to contain motif binding
struct MotifMap{
	map<string, GBMap> e;
};

// A structure to store prediction results
typedef struct{
	int TP;
	int FP;
	int TN;
	int FN;
} Pred;

// A structure to store prediction results as in Beer's format.
typedef struct{
        int label;
        double prob;
        string name;
} BPred;

// A structure to store the best solution.
struct BSolu{
	double s;
	vector<Constraint> cons;
	vector<CPTRow> cpt;
	set<int> mbnd;
	vector<MotifScore> mscor;
};

#endif


