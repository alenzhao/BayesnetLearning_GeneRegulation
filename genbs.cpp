#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include "bkgsub.h"

using namespace std;

int main()
{
	srand((unsigned)time(NULL));

	// Input node gene list.
	cout << "Input node gene list file: ";
	string node;
	cin >> node;
	vector<string> nlst;
	if(get1stcol(node, nlst) < 0)
		return 1;
	
	// Input bkg gene list.
	cout << "Input bkg gene list file: ";
	string bkg;
	cin >> bkg;
	vector<string> blst;
	if(get1stcol(bkg, blst) < 0)
		return 1;
	
	cout << "Input number of bootstrap samples: ";
	int bsn;
	cin >> bsn;	
	// Do 10 bootstrap sampling.
	for(int bs = 0; bs < bsn; bs++)
	{
		ostringstream strmo;
		strmo << node << bs+1;
		ofstream o(strmo.str().data());
		vector<string> rlst;
		for(size_t i = 0; i < nlst.size(); i++)
		{
			int idx = (int)floor(((double)rand()-1)/RAND_MAX*nlst.size());
			rlst.push_back(nlst[idx]);
		}

		for(size_t i = 0; i < rlst.size(); i++)
			o << rlst[i] << endl;
		o.close();
	}

	for(int bs = 0; bs < bsn; bs++)
	{
		ostringstream strmo;
		strmo << bkg << bs+1;
		ofstream o(strmo.str().data());
		vector<string> rlst;
		for(size_t i = 0; i < blst.size(); i++)
		{
			int idx = (int)floor(((double)rand()-1)/RAND_MAX*blst.size());
			rlst.push_back(blst[idx]);
		}

		for(size_t i = 0; i < rlst.size(); i++)
			o << rlst[i] << endl;
		o.close();
	}

	return 0;
}












