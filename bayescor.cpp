/*	bayescor.cpp

	1. Load a cluster's motif list, motifs' binding and gene lists into memory.
	2. Run Bayesian network to learn presense for each motif individually.
	3. Output the score for each motif.
*/


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "bayesub.h"
#include "globals.h"
#include "CmdLine.h"


int main(int argc, char* argv[])
{
	CCmdLine cmdLine;

	if(cmdLine.SplitLine(argc, argv) < 5)
	{
		cerr << "Usage: ./bayescor -m motif_list -n node_list -b bkg_list -f func_depth_folder -o output" << endl;
		cerr << "-i\tUse mutual information instead of Bayesian score" << endl;
		cerr << endl << "This calculate single motif's presence score on a cluster." << endl;
		cerr << "You need to run it before BBNet & GBNet." << endl;
		cerr << endl << "Contact: \"Li Shen\"<shen@ucsd.edu>" << endl;
		return 1;
	}

	string m, n, b, f, o;
	try
	{
		m = cmdLine.GetArgument("-m", 0);
		n = cmdLine.GetArgument("-n", 0);
		b = cmdLine.GetArgument("-b", 0);
		f = cmdLine.GetArgument("-f", 0);
		o = cmdLine.GetArgument("-o", 0);
	}
	catch(int)
	{
		cerr << "Wrong arguments!" << endl;
		return 1;
	}

	if(cmdLine.HasSwitch("-i"))
		itag = true;

	//itag = true;
	//string m = "../gbnet/data/Beer/motifs.list";
	//string n = "../gbnet/data/Beer/node.list";
	//string b = "../gbnet/data/Beer/bkg.list";
	//string f = "../gbnet/data/Beer/func";
	//string o = "../gbnet/data/Beer/scor_test.list";

	// Load motif list.
	vector<string> motiflst;
	if(loadmotif(motiflst, m) != 0)
	{
		cerr << "Load motif error!" << endl;
		return 1;
	}
	else
	{
#ifdef VERBOSE
		cout << "Load motif list completed!" << endl;
#endif
	}

	// Load gene list.
	vector<Case> tlst, blst;
	set<string> genset;	// Gene set for index.
	if(loadgene(tlst, blst, n, b) != 0)
	{
		cerr << "Load gene lists error!" << endl;
		return 1;
	}
	else
	{
#ifdef VERBOSE
		cout << "Load gene list completed!" << endl;
#endif
		for(size_t i = 0; i < tlst.size(); i++)
			genset.insert(tlst[i].name);
		for(size_t i = 0; i < blst.size(); i++)
			genset.insert(blst[i].name);
	}
	vector<Case> genlst;
	genlst.insert(genlst.end(), tlst.begin(), tlst.end());
	genlst.insert(genlst.end(), blst.begin(), blst.end());

	// Load motif binding information.
	if(loadbind(allbind, motiflst, genset, f) != 0)
	{
		cerr << "Load binding information error!" << endl;
		return 1;
	}
	else
	{
#ifdef VERBOSE
		cout << "Load binding information completed!" << endl;
#endif
	}

	mbnd.clear();
	mbnd.insert(0);	// Always only one motif in the list.
	vector<MotifScore> vscor;	// Store the final results for all motifs.
	// Calculate Bayesian score for each motif at each functional depth.
	for(size_t i = 0; i < motiflst.size(); i++)
	{
		mscor.clear();	// Clear all motifs' score and depth to restart.
		MotifScore scor;	// scor is used to initialize mscor and store the best.
		scor.name = motiflst[i];
		scor.score = 1.0;
		mscor.push_back(scor);

#ifdef VERBOSE
		cout << "Calculating score for motif " << motiflst[i] << "..." << endl;
#endif
		for(int j = 0; j < nfunc; j++)
		{
#ifdef VERBOSE
			cout << "Choosing functional depth " << func_depths[j] << "..." << endl;
#endif
			vector<CPTRow> cpt, ppt;

			Constraint pres = {"pres", 0, -1, -1};
			vector<Constraint> cons;
			cons.push_back(pres);

			mscor[0].depth = func_depths[j];
			constrcpt(cpt, ppt, genlst, cons);
			double s;
			if(!itag)
				s = score(1, cpt, ppt);
			else
				s = iscore(1, cpt);
			if(scor.score == 1 || s > scor.score)
			{
				scor.score = s;
				scor.depth = func_depths[j];
			}
		}
		vscor.push_back(scor);
	}

	// Output all motifs' scores and optimal depths.
	ofstream hScor(o.data());
	if(!hScor)
	{
		cerr << "Can't open " << o << endl;
		return 1;
	}
	outscor(hScor, vscor);
	hScor.close();

	return 0;
}


