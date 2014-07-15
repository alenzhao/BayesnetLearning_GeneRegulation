/*	newfunc.cpp

	Do pre-processing step for Bayesian network.
	Generate scores of weight matrices on each gene's sequence.
	Use the score files to generate function depth files.
*/


#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <assert.h>
#include "prepsub.h"
#include "CmdLine.h"

using namespace std;

const double fthrld = 0.01;

int main(int argc, char* argv[])
{
	// Read in parameters from command line.
	CCmdLine cmdLine;

	if(cmdLine.SplitLine(argc, argv) < 4)
	{
		cerr << "Usage: ./func -m motif_list -w pwm_folder -g genomic_sequence -f func_folder [-n normalization]" << endl;
		return 1;
	}

	string m, w, g, f;
	try
	{
		m = cmdLine.GetArgument("-m", 0);	// Motif list.
		w = cmdLine.GetArgument("-w", 0);	// PWM folder.
		g = cmdLine.GetArgument("-g", 0);	// Genomic sequence.
		f = cmdLine.GetArgument("-f", 0);	// Functional depth folder.
	}
	catch(int)
	{
		cerr << "Wrong arguments!" << endl;
		return 1;
	}

	string n = cmdLine.GetSafeArgument("-n", 0, "");	// Output motifs and corresponding normalization scores.

	// Read in motif list.
	vector<string> tf_list;
	if(get1stcol(m, tf_list) < 0)
		return 1;

	// File for storing all normalizing constants.
	ofstream hNorm;
	if(n != "")
	{
		hNorm.open(n.data());
		if(!hNorm)
		{
			cerr << "Can't open " << n << endl;
			return 1;
		}
	}
	
	// Load sequences.
	GenMap genmap;
	if(getseq(g, genmap) != 0)
		return 1;

	for(unsigned loop = 0; loop < tf_list.size(); loop++)
	{
		string tf = tf_list[loop];
		cout << "Processing TF: " << tf << endl;
		// Read weight matrix into memory.
		string fTF = w + "/" + tf + ".wm";
		double wm[MAXLEN][4];
		int mLen = loadpwm(fTF, wm);
		if(mLen < 0)
			return 1;

		// Read one gene at a time and calculate the scores of the weight matrix at both directions.
		double norm_const = MINVAL; // normalizing constant for this TF.
		for(GenMap::iterator g = genmap.begin(); g != genmap.end(); g++)
		{
			const string& gene = g->first;	// gene name.
			const string& seq = g->second.seq;	// sequence.
			size_t pLen = seq.length(); // promoter length.
			vector<double>& fscor = g->second.fscor;
			fscor.resize(pLen - mLen + 1);
			// Forward direction first.
			for(size_t i = 0; i < pLen - mLen + 1; i++)
			{
				double score = matscore(seq, i, wm, mLen);	// calculate the score at one position.
				fscor[i] = score;	// store the score in a vector.
				if(score > norm_const)
					norm_const = score;
			}
			// Reverse direction.
			string rev = copyrev(seq);
			vector<double>& rscor = g->second.rscor;
			rscor.resize(pLen - mLen + 1);
			for(size_t i = 0; i < pLen - mLen + 1; i++)
			{
				double score = matscore(rev, i, wm, mLen);
				rscor[pLen-mLen-i] = score;
				if(score > norm_const)
					norm_const = score;
			}
		}
		// Normalize both forward and reverse score matrices.
		for(GenMap::iterator i = genmap.begin(); i != genmap.end(); i++)
		{
			size_t pLen = i->second.seq.length();
			for(size_t j = 0; j < pLen - mLen + 1; j++)
			{
				i->second.fscor[j] /= norm_const;
				i->second.rscor[j] /= norm_const;
			}
		}
		if(n != "")
			hNorm << tf << '\t' << norm_const << endl;	// store normalizing constant for record.

		// Output func depth file.
		ostringstream sfunc;
		sfunc << f << "/" << tf << ".func";
		string ffunc = sfunc.str();
		ofstream hfunc(ffunc.data());
		if(!hfunc)
		{
			cerr << "Can't open " << ffunc << endl;
			return 1;
		}

		for(GenMap::iterator i = genmap.begin(); i != genmap.end(); i++)
		{
			hfunc << i->first;
			size_t pLen = i->second.seq.size();
			ostringstream sline;	// a line contain all func depth info for a gene.
			sline.precision(2);
			const vector<double>& fscor = i->second.fscor;
			const vector<double>& rscor = i->second.rscor;
			assert(fscor.size() == rscor.size());	// forward and reverse score vector should have same length.
			int count = 0;
			for(size_t j = 0; j < fscor.size(); j++)
			{
				if(fscor[j] >= fthrld)
				{
					sline << "\tF," << fscor[j] << ',' << pLen-mLen-j;
					count++;
				}
				if(rscor[j] >= fthrld)
				{
					sline << "\tR," << rscor[j] << ',' << pLen-mLen-j;
					count++;
				}
			}
			hfunc << "\t" << count << sline.str() << endl;
		}
	}

	return 0;
}










