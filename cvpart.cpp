#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "bkgsub.h"
#include "CmdLine.h"

int main(int argc, char* argv[])
{
	// Read parameters from console.
	CCmdLine cmdLine;

	if(cmdLine.SplitLine(argc, argv) < 3)
	{
		cerr << "Not enough arguments!" << endl;
		return 1;
	}

	string f, n, o;
	try
	{
		f = cmdLine.GetArgument("-f", 0);	// input folder.
		n = cmdLine.GetArgument("-n", 0);	// number of folds.
		o = cmdLine.GetArgument("-o", 0);	// output folder.
	}
	catch(int)
	{
		cerr << "Wrong arguments!" << endl;
		return 1;
	}

	string folder = f;
	int cvn = atoi(n.data());
	// Load and separate node gene list.
	string node = folder + "/node.list";
	vector<string> nlst;
	if(get1stcol(node, nlst) < 0)
		return 1;
	for(int i = 0; i < cvn; i++)
	{
		ostringstream strmPart;
		strmPart << o << "/node_part" << i+1 << ".list";
		string fPart = strmPart.str();
		ofstream hPart(fPart.data());
		if(!hPart)
		{
			cerr << "Can't open " << fPart << endl;
			return 1;
		}
		ostringstream strmTrn;
		strmTrn << o << "/node_trn" << i+1 << ".list";
		string fTrn = strmTrn.str();
		ofstream hTrn(fTrn.data());
		if(!hTrn)
		{
			cerr << "Can't open " << fTrn << endl;
			return 1;
		}
		for(int j = 0; j < nlst.size(); j++)
		{
			if((j-i)%cvn == 0)
				hPart << nlst[j] << endl;
			else
				hTrn << nlst[j] << endl;
		}
		hPart.close();
		hTrn.close();
	}
	
	// Load and separate bkg gene list.
	string bkg = folder + "/bkg.list";
	vector<string> blst;
	if(get1stcol(bkg, blst) < 0)
		return 1;
	for(int i = 0; i < cvn; i++)
	{
		ostringstream strmPart;
		strmPart << o << "/bkg_part" << i+1 << ".list";
		string fPart = strmPart.str();
		ofstream hPart(fPart.data());
		if(!hPart)
		{
			cerr << "Can't open " << fPart << endl;
			return 1;
		}
		ostringstream strmTrn;
		strmTrn << o << "/bkg_trn" << i+1 << ".list";
		string fTrn = strmTrn.str();
		ofstream hTrn(fTrn.data());
		if(!hTrn)
		{
			cerr << "Can't open " << fTrn << endl;
			return 1;
		}
		for(int j = 0; j < blst.size(); j++)
		{
			if((j-i)%cvn == 0)
				hPart << blst[j] << endl;
			else
				hTrn << blst[j] << endl;
		}
		hPart.close();
		hTrn.close();
	}

	return 1;
}








