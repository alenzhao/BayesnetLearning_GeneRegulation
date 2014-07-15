/*	corrbkg.cpp

	Find background nodes based on Pearson correlation.
*/


#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <string>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <algorithm>
#include "bkgsub.h"

using namespace std;

typedef struct{
	string name;
	double corr;
} NCorr;

bool cmp(NCorr a, NCorr b)
{
	return a.corr > b.corr;
}

int main()
{
	cout << "Input foreground cluster file: ";
	string node;
	cin >> node;
	
	cout << "Input folder for all background clusters: ";
	string folder;
	cin >> folder;
	
	// 158 for GNF Human; 122 for GNF Mouse.
	cout << "Number of columns: ";
	int cols;
	cin >> cols;
	
	cout << "Threshold of Pearson correlation: ";
	double thrld;
	cin >> thrld;
	
	cout << "Genome expression file: ";
	string fExpr;
	cin >> fExpr;
	
	// Load node name list.
	string fList = folder + "/BGNnh.txt";
	vector<string> nlst;
	if(get1stcol(fList, nlst) < 0)
		return 1;

	// Load the foreground cluster gene list.
	vector<string> nglst;
	if(get1stcol(node, nglst) < 0)
	{
		cerr << "Load cluster gene list error!" << endl;
		return 1;
	}

	// Load all gene expression profiles.
	map<string, vector<double> > exprmap;
	if(loadexpr(fExpr, exprmap, cols) != 0)
	{
		cerr << "Load expression profiles error!" << endl;
		return 1;
	}

	// Calculate the mean expression profile of the cluster.
	vector<double> mv = meanprof(nglst, exprmap);

	// Find background nodes. Output node list.
	cout << "Output background gene file: ";
	string fBG;
	cin >> fBG;
	ofstream hBG(fBG.data());
	if(!hBG)
	{
		cerr << "Can't open " << fBG << endl;
		return 1;
	}
	vector<NCorr> ncorr;
	for(size_t i = 0; i < nlst.size(); i++)
	{
		string name = nlst[i];
		if(name == node)
			continue;
		string file = folder + "/" + name + ".lst";
		vector<string> bglst;
		if(get1stcol(file, bglst) < 0)
		{
			cerr << "Load background node " << file << " error!" << endl;
			return 1;
		}
		vector<double> mb = meanprof(bglst, exprmap);
		double corr = Pearcorr0(mv, mb);
		NCorr one;
		one.name = name;
		one.corr = corr;
		ncorr.push_back(one);
		if(corr < thrld)
		{
			for(size_t j = 0; j < bglst.size(); j++)
			{
				//if(Pearcorr(mv, exprmap[bglst[j]]) < thrld)
					hBG << bglst[j] << endl;
			}
		}
	}
	hBG.close();

	cout << "Output background cluster file (rank in correlation): ";
	string fOut;
	cin >> fOut;
	ofstream hOut(fOut.data());
	if(!hOut)
	{
		cerr << "Can't open " << fOut << endl;
		return 1;
	}
	sort(ncorr.begin(), ncorr.end(), cmp);
	for(size_t i = 0; i < ncorr.size(); i++)
		hOut << ncorr[i].name << "\t" << ncorr[i].corr << endl;
	hOut.close();


	return 0;
}


