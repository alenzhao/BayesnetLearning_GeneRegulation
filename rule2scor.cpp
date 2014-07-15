/*	rule2scor.cpp

	Use a set of rules to calculate the Bayesian score of a set of genes.
*/


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include "mysub.h"
#include "CmdLine.h"


// Functional depths.
const double func_depths[] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 
						0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95};
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
int motifcand = 5;	// maximum number of motifs to be included.
double logK = 1.5;		// parameter for network structure prior.
int prior = 0;	// flag for setting prior counts for motifs that should be included into Bayesian networks.
map<string, bool> primo;	// map to store motifs that should be added to Bayesian network apriori.
int pricnt = 20;	// prior counts for preferred motifs.
int bootstrap = -1;	// Index for bootstrap sample.
int bsn = 10;	// number of bootstrap samples.
bool usebs = true;	// Use bootstrapping or not.

// Load constraints and CPT.
int loadcons(const string& c, vector<Constraint>& cons)
{
	ifstream h(c.data());
	if(!h)
	{
		cerr << "Can't open " << c << endl;
		return 1;
	}
	while(h.good())
	{
		string strLn;
		getline(h, strLn);
		if(strLn == "")
			continue;
		istringstream strmLn(strLn);
		Constraint con;
		strmLn >> con.desc >> con.motif0 >> con.motif1 >> con.para;
		cons.push_back(con);
	}

	return 0;
}

// Extract binding information.
void extrbnd(const vector<MotifScore>& mscor, vector<MBinding>& mbnd, MotifMap& allbind)
{
	for(int i = 0; i < mscor.size(); i++)
	{
		string motif = mscor[i].name;
		double depth = mscor[i].depth;
		int didx = depidx(depth);
		MBinding onebind;
		onebind.name = motif;
		onebind.depth = depth;
		onebind.bindinfo = allbind.e[motif].e[didx];
		mbnd.push_back(onebind);
	}
}

int main(int argc, char* argv[])
{
	// Read parameters from console.
	CCmdLine cmdLine;

	if(cmdLine.SplitLine(argc, argv) < 6)
	{
		cerr << "Not enough arguments!" << endl;
		return 1;
	}

	string m, c, n, b, f, k;
	double logK;
	try
	{
		m = cmdLine.GetArgument("-m", 0);	// motif file.
		c = cmdLine.GetArgument("-c", 0);	// constraints and CPT file.
		n = cmdLine.GetArgument("-n", 0);	// node gene list file.
		b = cmdLine.GetArgument("-b", 0);	// bkg gene list file.
		f = cmdLine.GetArgument("-f", 0);	// func depth folder.
		k = cmdLine.GetArgument("-k", 0);	// logK value.
		logK = atof(k.data());
	}
	catch(int)
	{
		cerr << "Wrong arguments!" << endl;
		return 1;
	}
	
	// Load motif Bayesian score file.
	vector<MotifScore> mscor;
	if(loadscor(mscor, m) != 0)
	{
		cerr << "Load motif list error!" << endl;
		return 1;
	}

	// Load gene list.
	vector<Case> tlst, blst, genlst;
	map<string, bool> genmap;
	if(loadgene(tlst, blst, n, b) != 0)
	{
		cerr << "Load gene lists error!" << endl;
		return 1;
	}
	else
	{
		genlst.insert(genlst.end(), tlst.begin(), tlst.end());
		genlst.insert(genlst.end(), blst.begin(), blst.end());
		cout << "Load gene list completed!" << endl;
		// All training and testing gene names are put into genmap.
		for(size_t i = 0; i < genlst.size(); i++)
			genmap[genlst[i].name] = true;
	}

	// Load motif binding information of genes in genmap.
	MotifMap allbind;
	if(loadbind(allbind, mscor, genmap, f) != 0)
	{
		cerr << "Load binding information error!" << endl;
		return 1;
	}
	else
		cout << "Load binding information completed!" << endl;
			

	vector<Constraint> cons;
	loadcons(c, cons);	// Load constraints.
	vector<MBinding> mbnd;
	extrbnd(mscor, mbnd, allbind);	// Extract binding info.
	vector<CPTRow> cpt, ppt;
	constrcpt(cpt, ppt, genlst, cons, mbnd);	// Construct CPT.

	double s = score(cons.size(), cpt, ppt);	// Calculate Bayesian score.
	cout << "Bayesian score = " << s << endl;

	return 0;
}



