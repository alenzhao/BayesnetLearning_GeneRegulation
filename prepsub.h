#include <sstream>
#include <string>
#include <assert.h>
#include <math.h>

using namespace std;


double normal_01_cdf ( double x );
double normal_cdf ( double x, double a, double b );


//******************************************************************************

double normal_01_cdf ( double x )

//******************************************************************************
//
//  Purpose:
//
//    NORMAL_01_CDF evaluates the Normal 01 CDF.
//
//  Modified:
//
//    10 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference: 
//
//    A G Adams,
//    Areas Under the Normal Curve,
//    Algorithm 39, 
//    Computer j., 
//    Volume 12, pages 197-198, 1969.
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Output, double CDF, the value of the CDF.
//
{
	double a1 = 0.398942280444;
	double a2 = 0.399903438504;
	double a3 = 5.75885480458;
	double a4 = 29.8213557808;
	double a5 = 2.62433121679;
	double a6 = 48.6959930692;
	double a7 = 5.92885724438;
	double b0 = 0.398942280385;
	double b1 = 3.8052E-08;
	double b2 = 1.00000615302;
	double b3 = 3.98064794E-04;
	double b4 = 1.98615381364;
	double b5 = 0.151679116635;
	double b6 = 5.29330324926;
	double b7 = 4.8385912808;
	double b8 = 15.1508972451;
	double b9 = 0.742380924027;
	double b10 = 30.789933034;
	double b11 = 3.99019417011;
	double cdf;
	double q;
	double y;
	//
	//  |X| <= 1.28.
	//
	if ( fabs ( x ) <= 1.28 )
	{
		y = 0.5 * x * x;

		q = 0.5 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 
			+ a6 / ( y + a7 ) ) ) );
		//
		//  1.28 < |X| <= 12.7
		//
	}
	else if ( fabs ( x ) <= 12.7 )
	{
		y = 0.5 * x * x;

		q = exp ( - y ) * b0 / ( fabs ( x ) - b1 
			+ b2 / ( fabs ( x ) + b3 
			+ b4 / ( fabs ( x ) - b5 
			+ b6 / ( fabs ( x ) + b7 
			- b8 / ( fabs ( x ) + b9 
			+ b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
		//
		//  12.7 < |X|
		//
	}
	else
	{
		q = 0.0;
	}
	//
	//  Take account of negative X.
	//
	if ( x < 0.0 )
	{
		cdf = q;
	}
	else
	{
		cdf = 1.0 - q;
	}

	return cdf;
}
//******************************************************************************

double normal_cdf ( double x, double a, double b )

//******************************************************************************
//
//  Purpose:
//
//    NORMAL_CDF evaluates the Normal CDF.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double CDF, the value of the CDF.
//
{
	double cdf;
	double y;

	y = ( x - a ) / b;

	cdf = normal_01_cdf ( y );

	return cdf;
}


// ************
// Start non-normal routines.
// ************

#define CUTOFF 1e-3
#define NDCUT 5
#define NACUT 5
#define RESAMP 1000
#define MAXLEN	40
#define MINVAL -9999.0

typedef struct{
	char orien;
	int loc;
	long double score;
} BindInfo;

typedef struct {
	char orien1, orien2;
	int loc1, loc2;
} BInfo2;

typedef struct {
	double mean;
	double std;
} NormPara;

typedef struct{
	string seq;	// promoter sequence.
	vector<double> fscor;	// forward score vector.
	vector<double> rscor;	// reverse score vector.
} GenInfo;	// Score information for one gene.
typedef map<string, GenInfo> GenMap;	// Map to store all genes' score information.

unsigned dnaidx(char base); // Give an index to DNA base.
char basepair(char base); // base-pairing function.

// Transform a string to upper case.
string& str2upper(string& str)
{
	for(unsigned int i = 0; i < str.length(); i++)
		str[i] = toupper(str[i]);

	return str;
}

// Calculate the mean of a vector.
template<class T>
double mean(const vector<T>& v)
{
	double m = 0.0;
	for(unsigned i = 0; i < v.size(); i++)
		m += v[i];
	m /= v.size();

	return m;
}

// Calculate the standard deviation of binding location.
template<class T>
double stnd(const vector<T>& v, double m)
{
	double ss = 0.0;
	for(unsigned i = 0; i < v.size(); i++)
		ss += (v[i] - m)*(v[i] - m);
	ss /= (v.size()-1);

	return sqrt(ss);
}

// Read binding information of one TF and produce a hash map.
int getfbind(const string& tf, const string& spec, map<string, BindInfo>& bMap)
{
	string fTF = "results/" + spec + "_Bind/" + tf + ".bind";
	ifstream hTF(fTF.data());
	if(!hTF)
	{
		cerr << "Can't open " << fTF << endl;
		return 1;
	}

	while(hTF.good())
	{
		BindInfo bInfo;
		string gene;
		hTF >> gene;
		if(gene == "")
			continue;
		str2upper(gene);
		hTF >> bInfo.orien;
		hTF >> bInfo.loc;
		hTF >> bInfo.score;
		hTF.ignore(1, '\n');

		bMap[gene] = bInfo;
	}
	hTF.close();

	return 0;
}

// Get the binding information of the two TFs on one cluster.
void getnbind(const vector<string>& node, map<string, BindInfo>& bMap1, 
			  map<string, BindInfo>& bMap2, map<string, BInfo2>& nMap)
{
	for(unsigned i = 0; i < node.size(); i++)
	{
		BindInfo bind1 = bMap1[node[i]];
		BindInfo bind2 = bMap2[node[i]];
		if(bind1.score > CUTOFF || bind2.score > CUTOFF)
			continue;
		BInfo2 bInfo2;
		bInfo2.orien1 = bind1.orien;
		bInfo2.loc1 = bind1.loc;
		bInfo2.orien2 = bind2.orien;
		bInfo2.loc2 = bind2.loc;
		nMap[node[i]] = bInfo2;
	}
}

// Get the node gene list.
int getnglst(const string& node, vector<string>& nglst)
{
	string fNode = "data/NODE" + node + ".lst";
	ifstream hNode(fNode.data());
	if(!hNode)
	{
		cerr << "Can't open " << fNode << endl;
		return 1;
	}

	while(hNode.good())
	{
		string gene;
		hNode >> gene;
		hNode.ignore(1, '\n');
		if(gene == "")
			continue;
		nglst.push_back(gene);
	}
	hNode.close(); // node gene list read finished!

	return 0;
}

const int dcutoff[NDCUT] = {100, 200, 300, 400, 500};
const int atgcutoff[NACUT] = {100, 200, 300, 400, 500};

typedef struct{
	string vname;
	unsigned len;
	ifstream hf;
	ifstream hr;
	double mean;
	double stnd;
} TFInfo;

typedef struct{
	double mean;
	double stnd;
} StatInfo;

// Read in one line of socres.
string readscor(ifstream& h, vector<double>& scores, const unsigned mLen, const string& cgen = "")
{
	string strLn;
	getline(h, strLn);
	if(strLn == "")
		return "";
	istringstream lnStream(strLn);
	string gene;
	unsigned pLen;
	lnStream >> gene >> pLen;
	str2upper(gene);
	if(cgen == gene || cgen == "")
	{
		scores.resize(pLen - mLen);
		for(unsigned i = 0; i < pLen - mLen; i++)
			lnStream >> scores[i];
	}

	return gene;
}

// Read in hash table for TF versions and their lengths.
int readtflen(const string& f, map<string, unsigned>& m)
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return 1;
	}

	while(h.good())
	{
		string tf;
		unsigned len;
		h >> tf >> len;
		if(tf == "")
			continue;

		m[tf] = len;
	}
	h.close();

	return 0;
}

// Read in hash table for TF versions and their means and standard deviations.
int readtfstat(const string& f, map<string, StatInfo>& m)
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return 1;
	}
	while(h.good())
	{
		string tf;
		StatInfo tfInfo;
		h >> tf >> tfInfo.mean >> tfInfo.stnd;
		h.ignore(1, '\n');
		if(tf == "")
			continue;
		m[tf] = tfInfo;
	}
	h.close();

	return 0;
}

// Data structure to store all TF versions' names, lengths and file handles.
int getfinfo(TFInfo infoVec[], const string& tfName, const string& spec, map<string, vector<string> >& verMap, 
			 map<string, unsigned>& TFLen, map<string, StatInfo>& TFStat)
{
	vector<string>& vers = verMap[tfName];
	for(unsigned i = 0; i < vers.size(); i++)
	{
		infoVec[i].vname = vers.at(i);
		string fFor = "scores/" + spec + "/" + vers[i] + ".forward.scores";
		infoVec[i].hf.open(fFor.data());
		string fRev = "scores/" + spec + "/" + vers[i] + ".reverse.scores";
		infoVec[i].hr.open(fRev.data());
		infoVec[i].len = TFLen[vers[i]];
		infoVec[i].mean = TFStat[vers[i]].mean;
		infoVec[i].stnd = TFStat[vers[i]].stnd;
	}

	return (int)vers.size();
}

// Convert the binding scores to binding specificity using negative P-values.
void scor2bind(vector<double>& scor, double m, double s)
{
	for(unsigned i = 0; i < scor.size(); i++)
	{
		double pval = 1 - normal_cdf(scor.at(i), m, s);
		scor.at(i) = -log10(pval);
	}
}

void combimax(vector<double>& bind, const vector<double>& scor)
{
	if(bind.empty())
		bind = scor;
	else
	{
		if(scor.size() > bind.size())
		{
			bind.insert(bind.begin(), scor.begin(), scor.begin()+scor.size()-bind.size());
			assert(scor.size() == bind.size());
		}
		int bIdx = (int)bind.size()-1;
		int sIdx = (int)scor.size()-1;
		for(unsigned i = 0; i < scor.size(); i++)
		{
			if(bind.at(bIdx) < scor.at(sIdx))
				bind.at(bIdx) = scor.at(sIdx);
			bIdx--;
			sIdx--;
		}
	}
}

// Get the binding specificity of a TF on a gene.
void getfgenbind(vector<double>& fBindSpec, vector<double>& rBindSpec, TFInfo infoVec[], 
				 const vector<string>& vers, const string& gene)
{
	bool tag = false;	// Tag for determining whether the gene has been read.
	while(infoVec[0].hf.good())
	{
		for(unsigned i = 0; i < vers.size(); i++)
		{
			vector<double> fscor;
			string fgen = readscor(infoVec[i].hf, fscor, infoVec[i].len, gene);
			if(fgen == "")
				break;
			if(fgen != gene)
			{
				readscor(infoVec[i].hr, fscor, infoVec[i].len, gene);	// read a NULL line.
				continue;
			}
			tag = true;
			scor2bind(fscor, infoVec[i].mean, infoVec[i].stnd);
			combimax(fBindSpec, fscor);

			vector<double> rscor;
			string rgen = readscor(infoVec[i].hr, rscor, infoVec[i].len, gene);
			assert(fgen == rgen);
			scor2bind(rscor, infoVec[i].mean, infoVec[i].stnd);
			combimax(rBindSpec, rscor);
		}
		if(tag)
			break;
	}
}

// Given a file, read in the first column as a vector of strings.
int get1stcol(const string& f, vector<string>& list)
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return -1;
	}

	while(h.good())
	{
		string strLn;
		getline(h, strLn);
		if(strLn == "")
			continue;
		istringstream strmLn(strLn);
		string elem;
		strmLn >> elem;
		list.push_back(elem);
	}
	h.close();

	return (int)list.size();
}

char basepair(char base) // base-pairing function.
{
	switch(base)
	{
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	default:
		return 'N';
	}
}

// Copy reversed sequence into another string.
string copyrev(const string& seq)
{
	string rev = seq;
	for(unsigned int i = 0; i < seq.length(); i++)
		rev[i] = basepair(seq[seq.length()-1-i]);

	return rev;
}

// Calculate the score of a weight matrix at a position on one sequence.
double matscore(const string& seq, const size_t pos, const double wm[][4], const int len)
{
	double score = 1.0;
	for(size_t i = pos; i < pos + len; i++)
	{
		if(seq[i] == 'N')
			return MINVAL;
		score *= wm[i-pos][dnaidx(seq[i])];
	}

	return score;
}

unsigned dnaidx(char base) // Give an index to DNA base.
{
	switch(base)
	{
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
		return 99;
	}
}

// Load promoter sequences for all genes.
int getseq(const string& f, GenMap& m)
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << "!" << endl;
		return 1;
	}
	while(h.good())
	{
		string strLn;
		getline(h, strLn);
		if(strLn == "")
			continue;
		size_t tab = strLn.find('\t');
		string gene = strLn.substr(0, tab);
		GenInfo info;
		info.seq = strLn.substr(tab+1);
		m[gene] = info;
	}
	return 0;
}

// Load a PWM from file. Return length of positions.
int loadpwm(const string& f, double wm[MAXLEN][4])
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return -1;
	}

	int len = 0, pos;
	string headline;
	getline(h, headline);	// Get rid of the header line.
	while(h.good())
	{
		string strLn;
		getline(h, strLn);
		if(strLn == "")
			continue;

		istringstream sline(strLn);
		double A, C, G, T;
		sline >> pos >> A >> C >> G >> T;
		if(pos != len)
			cerr << "Error in position index of PWM file! Ignore..." << endl;

		wm[len][0] = A;
		wm[len][1] = C;
		wm[len][2] = G;
		wm[len][3] = T;
		len++;
	}
	return len;
}


