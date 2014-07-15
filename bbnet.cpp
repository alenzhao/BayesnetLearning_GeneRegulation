/*	bbnet.cpp

	Main routine of Bayesian network.
	1. Load data into memory.
	2. Run Bayesian network and infer the rules.
	3. Output the results.
*/


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include "bayesub.h"
#include "globals.h"
#include "CmdLine.h"


int main(int argc, char* argv[])
{
	// Read parameters from console.
	CCmdLine cmdLine;

	if(cmdLine.SplitLine(argc, argv) < 5)
	{
		cerr << "Usage: ./bbnet -s score_file -n node -b bkg -f func_depth -o output" << endl;
		cerr << endl << "Additional parameter:" << endl;
		cerr << "-k\tPenalty parameter(logK, Default = 5.0)" << endl;
		cerr << "-c\tnumber of candidate motifs (Default=50)" << endl;
		cerr << "-d\tpositive negative (for prediction using BN)" << endl;
		cerr << "-l\toutput of all training samples' information." << endl;
		cerr << "-t\ttranslational(transcriptional) start sites.(Default = right end)" << endl;
		cerr << "-rb\tbit-string to determine which rules to include.(Default = 111110)" << endl;
		cerr << "-i\tUse mutual information instead of Bayesian score" << endl;
		cerr << endl << "Contact: \"Li Shen\"<shen@ucsd.edu>" << endl;
		return 1;
	}

	string s, n, b, f, o;
	try
	{
		s = cmdLine.GetArgument("-s", 0);	// score file.
		n = cmdLine.GetArgument("-n", 0);	// node gene list file.
		b = cmdLine.GetArgument("-b", 0);	// bkg gene list file.
		f = cmdLine.GetArgument("-f", 0);	// func depth folder.
		o = cmdLine.GetArgument("-o", 0);	// results output file.
	}
	catch(int)
	{
		cerr << "Wrong arguments!" << endl;
		return 1;
	}

	if(cmdLine.HasSwitch("-i"))
		itag = true;

	//itag = true;
	//string s = "../gbnet/data/Beer/scor_test.list";
	//string n = "../gbnet/data/Beer/node.list";
	//string b = "../gbnet/data/Beer/bkg.list";
	//string f = "../gbnet/data/Beer/func";
	//string k = "0.015";
	//logK = atof(k.data());
	//string o = "../gbnet/data/Beer/bb_res_test2.txt";

	string k;	// Penalty parameter; logK value.
	if(!itag)
		k = cmdLine.GetSafeArgument("-k", 0, "5.0");
	else
		k = cmdLine.GetSafeArgument("-k", 0, "0.015");
	logK = atof(k.data());

	string c = cmdLine.GetSafeArgument("-c", 0, "50");	// number of candidate motifs. default = 50.
	motifcand = atoi(c.data());

	// Use prior counts for some motifs if specified.
	string p = cmdLine.GetSafeArgument("-p", 0, "0");	// prior counts.
	pricnt = atoi(p.data());
	if(pricnt > 0)	// read preferred motifs list from file.
	{
		prior = 1;
		string fPrim = cmdLine.GetSafeArgument("-p", 1, "primot.txt");
		vector<string> primv;
		if(get1stcol(fPrim, primv) < 0)
			return 1;
		for(size_t i = 0; i < primv.size(); i++)
			primo.insert(primv[i]);
	}

	// File names for positive, negative and left-out testing lists.
	string pos = cmdLine.GetSafeArgument("-d", 0, "");	// positive testing cases.
	string neg = cmdLine.GetSafeArgument("-d", 1, "");	// negative testing cases.
	string res = cmdLine.GetSafeArgument("-d", 2, "");	// left-out testing cases.
	vector<string> plst, nlst, rlst;	// positive, negative and left-out lists.
	if(pos != "" && neg != "")
	{
		if(get1stcol(pos, plst) < 0)
			return 1;
		if(get1stcol(neg, nlst) < 0)
			return 1;
	}
	if(res != "")
	{
		if(get1stcol(res, rlst) < 0)
			return 1;
	}

	// File for output of all training samples' information.
	string finfo = cmdLine.GetSafeArgument("-l", 0, "");

	// File to store all genes' translational/transcriptional start sites.
	string ftss = cmdLine.GetSafeArgument("-t", 0, "");
	if(ftss != "")
		loadtss(ftss, mtss);

	// A bit-string to determine which rules to include.
	rb = cmdLine.GetSafeArgument("-rb", 0, "111110");

	string bp = cmdLine.GetSafeArgument("-bp", 0, "");	// Output each gene's probability like in Beer's prediction.
	
	// Load motif Bayesian score file.
	if(loadscor(mscor, s) != 0)
	{
		cerr << "Load motif scores eror!" << endl;
		return 1;
	}
	else
	{
#ifdef VERBOSE
		cout << "Display candidate motifs that are loaded:" << endl;
		dispscor(mscor);
#endif
	}
	vector<MotifScore> oscor = mscor;	// Save an original copy of motif scores.

	// Load gene list.
	vector<Case> tlst, blst, genlst;
	set<string> genset;
	if(loadgene(tlst, blst, n, b) != 0)
	{
		cerr << "Load gene lists error!" << endl;
		return 1;
	}
	else
	{
		genlst.insert(genlst.end(), tlst.begin(), tlst.end());
		genlst.insert(genlst.end(), blst.begin(), blst.end());
#ifdef VERBOSE
		cout << "Load gene list completed!" << endl;
#endif
		// All training and testing gene names are put into genmap.
		for(size_t i = 0; i < genlst.size(); i++)
			genset.insert(genlst[i].name);
		for(size_t i = 0; i < plst.size(); i++)
			genset.insert(plst[i]);
		for(size_t i = 0; i < nlst.size(); i++)
			genset.insert(nlst[i]);
		for(size_t i = 0; i < rlst.size(); i++)
			genset.insert(rlst[i]);
	}

	// Load motif binding information of genes in genmap.
	if(loadbind(allbind, mscor, genset, f) != 0)
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

	// File for output.
	ofstream hOut(o.data());
	if(!hOut)
	{
		cerr << "Can't open " << o << endl;
		return 1;
	}
	hOut << "Number of genes in category 1: " << tlst.size() << endl;
	hOut << "Number of genes in category 0: " << blst.size() << endl << endl;

#ifdef VERBOSE
	cout << endl << "Running on original data." << endl;
#endif
	vector<Constraint> cons;
	vector<CPTRow> cpt;
	clock_t start = clock();
	double scor = bbnet(cons, cpt, genlst);
	clock_t finish = clock();
	if(outbayes(hOut, scor, cons, cpt, oscor, tlst.size(), blst.size()) != 0)
	{
		cerr << "Output Bayesian network results error!" << endl;
		return 1;
	}
	if(finfo != "")
	{
		if(outgene(finfo, tlst, blst, cons) != 0)
			cerr << "Output training samples' information error!" << endl;
		return 1;
	}
	if(pos != "" && neg != "")
	{
		Pred d = predict(cons, cpt, plst, nlst);
		outpred(hOut, d, n, b, pos, neg);
		if(bp != "")	// output each gene's probability being in this cluster if output file is specified.
		{
			ofstream hbp(bp.data());
			if(!hbp)
			{
				cerr << "Can't open " << bp << endl;
				return 1;
			}
			vector<BPred> trnbp = predict(cons, cpt, genlst, 0);	// probabilities for training genes.
			outpred(hbp, trnbp);
			vector<string> tstlst;	// probabilities for testing genes.
			tstlst.insert(tstlst.end(), plst.begin(), plst.end());	// positive testings.
			tstlst.insert(tstlst.end(), nlst.begin(), nlst.end());	// negative testings.
			vector<BPred> tstbp = predict(cons, cpt, tstlst, 1);
			outpred(hbp, tstbp);
			if(res != "")	// probabilities for left-out genes if the left-out file is specified.
			{
				vector<BPred> lefbp = predict(cons, cpt, rlst, -1);
				outpred(hbp, lefbp);
			}
			hbp.close();
		}
	}
	hOut << endl << "Bayesian network occupied CPU " << (double)(finish-start)/CLOCKS_PER_SEC << " seconds." << endl;
	hOut.close();

	return 0;
}



