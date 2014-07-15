/*	bayesub.h

	Declarations of all mathematical and I/O subroutines.
*/

#ifndef BAYESUB_H
#define BAYESUB_H

#include <map>
#include <vector>
#include <iostream>
#include "badefs.h"

using namespace std;


// ************ All subroutines start here *************

// Learn Bayesian network - BBNet.
double bbnet(vector<Constraint>& cons, vector<CPTRow>& cpt, const vector<Case>& genlst);

// Learn Bayesian network - GBNet.
double gbnet(vector<Constraint>& cons, vector<CPTRow>& cpt, const vector<Case>& genlst);

// According to a set of constraints, classify a gene into a category. 
// Different combinations of the constraints are described in the bits of an integer.
int classification(const string& gene, const vector<Constraint>& cons);

// Test whether a gene satisfies one constraint.
int test(const Constraint& c, const VGB& m0, double d0, const VGB& m1, double d1, int tss = 0);

// Calculate Bayesian score given CPT and priors.
double score(int np, const vector<CPTRow>& cpt, const vector<CPTRow>& ppt);

// Calculate Normalized Mutual Information given CPT.
double iscore(int np, const vector<CPTRow>& cpt);

// Load motif scores from file and sort them in descending order.
int loadscor(vector<MotifScore>& mscor, const string& s);

// Load gene list from file.
int loadgene(vector<Case>& tlst, vector<Case>& blst, const string& n, const string& b);

// Load motif list from file.
int loadmotif(vector<string>& motiflst, const string& f);

// Load one motif's binding.
int loadone(GBMap& onebind, const string& motif, const set<string>& genset, const string& folder);

// Load all motifs' binding.
int loadbind(MotifMap& allbind, const vector<string>& motiflst, const set<string>& genset, const string& folder);
int loadbind(MotifMap& allbind, const vector<MotifScore>& mscor, const set<string>& genset, const string& folder); // Overloading.

// Comparing routine for sorting motif scores.
bool cmp(MotifScore s0, MotifScore s1);

// Display the candidate motifs and their scores.
void dispscor(const vector<MotifScore>& mscor);

// Construct conditional probability table given gene list, constraints and motif binding.
void constrcpt(vector<CPTRow>& cpt, vector<CPTRow>& ppt, const vector<Case>& genlst, const vector<Constraint>& cons);

// Extract binding of a site from a string.
GBinding extrbnd(const string& s);

// Transform a string to upper case.
string& str2upper(string& str);

// Calculate the log Gamma value given a integer.
double logamma(int x);

// Add one new constraint into Bayesian network and update everything if necessary.
double addcons(vector<Constraint>& cons, vector<CPTRow>& cpt, Constraint c, double s, 
			   const vector<Case>& genlst, const int paraset[], int npara, bool jump = false);

// Delete constraint to improve score.
double delcons(vector<Constraint>& cons, vector<CPTRow>& cpt, double s, const vector<Case>& genlst);

// Format and output the results of Bayesian network on a cluster.
int outbayes(ofstream& hOut, double s, const vector<Constraint>& cons, const vector<CPTRow>& cpt, const vector<MotifScore>& vms, size_t node, size_t bkg);

// Output one constraint using file handle.
void outcons(ofstream& h, const Constraint& c, const vector<MotifScore>& mscor);

// Format a non-negative integer to a string using binary representation.
string fmtbinary(int n, size_t t);

// Update functional depth of one motif to improve score.
double updepth(int mi, vector<Constraint>& cons, vector<CPTRow>& cpt, double s, const vector<Case>& genlst, bool jump = false);

// Check whether a constraint has already been added.
bool chkcons(const vector<Constraint>& cons, const string& desc, int motif0, int motif1 = -1);

// Add a presence node into Bayesian network.
double addpres(int mi, vector<Constraint>& cons, vector<CPTRow>& cpt, double s, const vector<Case>& genlst, bool jump = false);

// Output all motif scores and optimal functional depths to file.
void outscor(ofstream& h, const vector<MotifScore>& mscor);

// Find the index in depth array according to the depth.
int depidx(double depth);

// Add prior information into CPT.
void setprior(vector<CPTRow>& ppt, const vector<Constraint>& cons);

// Initialize CPT.
void initcpt(vector<CPTRow>& cpt, size_t ns, int val = 0);

// The number of genes that satisfy each constraint from CPT.
vector<CPTRow> ebitcpt(const vector<CPTRow>& cpt, size_t nc);

// Take a bootstrap sample for a vector of objects.
vector<Case> bsamp(const vector<Case>& t, const vector<Case>& b);

// Calculate number of cases in each category.
int calcnum(const vector<Case>& v, int c);

// Check whether a character is a decimal number.
bool isnum(char c);

// Given a file, read in the first column as a vector of strings.
int get1stcol(const string& f, vector<string>& list);

// Using rules, CPT and binding to predict cases, return prediction results.
Pred predict(const vector<Constraint>& cons, const vector<CPTRow>& cpt, const vector<string>& plst, const vector<string>& nlst);
// Predict each gene's probability of being in this cluster and output the list as in Beer's prediction.
vector<BPred> predict(const vector<Constraint>& cons, const vector<CPTRow>& cpt, const vector<string>& genlst, int label);
// Operator overload for different parameter.
vector<BPred> predict(const vector<Constraint>& cons, const vector<CPTRow>& cpt, const vector<Case>& genlst, int label);

// Output prediction results.
int outpred(ofstream& h, Pred d, const string& node, const string& bkg, const string& pos, const string& neg);
// Output each gene's probability.
int outpred(ofstream& h, const vector<BPred>& bp);

// Record the best solution.
inline void bestsolu(double s, const vector<Constraint>& cons, const vector<CPTRow>& cpt, const set<int> mbnd, const vector<MotifScore>& mscor);

// Restart SA with best solution if bad condition happens.
double restart(double s, vector<Constraint>& cons, vector<CPTRow>& cpt, int& iter);


// Calculate P-value based on Fisher's exact test.
double fpval(double nm, double nn, double bm, double bn);

// Output all genes' TF binding site information.
int outgene(const string& f, const vector<Case>& n, const vector<Case>& b, const vector<Constraint>& cons);

// Output one gene's TF binding site information.
void outbind(ofstream& h, const string& gene, const vector<Constraint>& cons);

// The binding sites that satisfy one constraint.
string binds(const Constraint& c, const VGB& m0, double d0, const VGB& m1, double d1);

// Format a string describing a binding site.
string fmtsite(const GBinding& gb, bool parth = true);

// Load translationl/transcriptional start sites.
int loadtss(const string& f, map<string, int>& m);

// Logrithm base 2 of value x.
inline double log2(double x);

#endif


