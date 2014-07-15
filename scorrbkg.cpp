/*	scorrbkg.cpp

	Select single genes as background against a cluster.
*/


#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "bkgsub.h"

using namespace std;

typedef map<string, vector<double> > EXPRMAP;

typedef struct{
	string name;
	double corr;
} GCORR;

bool cmp(GCORR a, GCORR b)
{
	return a.corr < b.corr;
}

int main()
{
	// Load the cluster gene list.
	cout << "Input cluster gene list file: ";
	string node;
	cin >> node;
	vector<string> nglst;
	if(get1stcol(node, nglst) < 0)
		return 1;

	// Load all gene expression profiles.
	cout << "Input genome expression file: ";
	string fExpr;
	cin >> fExpr;
	cout << "Number of columns: ";
	int cols;
	cin >> cols;
	EXPRMAP exprmap;
	if(loadexpr(fExpr, exprmap, cols) != 0)
	{
		cerr << "Load expression profiles error!" << endl;
		return 1;
	}

	// Calculate the mean expression profile of the cluster.
	vector<double> mv = meanprof(nglst, exprmap);

	// Calculate all genes' correlation to the mean of the cluster excluding YY1's targets.
	cout << "Input YY1's target gene list file: ";
	string fTar;
	cin >> fTar;
	vector<string> tarlst;
	if(get1stcol(fTar, tarlst) < 0)
		return 1;
	map<string, bool> tarmap;
	for(size_t i = 0; i < tarlst.size(); i++)
		tarmap[tarlst[i]] = true;

	vector<GCORR> bkglst;
	for(EXPRMAP::const_iterator i = exprmap.begin(); i != exprmap.end(); i++)
	{
		if(tarmap.find(i->first) != tarmap.end())
			continue;
		GCORR gcorr;
		gcorr.name = i->first;
		gcorr.corr = Pearcorr0(i->second, mv);
		bkglst.push_back(gcorr);
	}
	sort(bkglst.begin(), bkglst.end(), cmp);

	// Output the background gene list to a file.
	cout << "Output background gene file: ";
	string fOut;
	cin >> fOut;
	ofstream hOut(fOut.data());
	if(!hOut)
	{
		cerr << "Can't open " << fOut << endl;
		return 1;
	}
	for(size_t i = 0; i < bkglst.size(); i++)
		hOut << bkglst[i].name << "\t" << bkglst[i].corr << endl;
	hOut.close();
	
	return 0;
}


