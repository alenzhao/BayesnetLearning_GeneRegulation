#include <assert.h>
#include <math.h>

using namespace std;

#define MISSING -99999.0

int getnglst(const string& node, vector<string>& nglst);
int getseq(const string& f, map<string, string>& m);
string& str2upper(string& str);
double Pearcorr(const vector<double>& x, const vector<double>& y);
double Pearcorr0(const vector<double>& x, const vector<double>& y);
double meancorr(const vector<string>& gset, map<string, vector<double> >& exprmap, bool uncen = false);
int loadexpr(const string& f, map<string, vector<double> >& exprMap, size_t cols);
template<class T>
void str2vec(const string& str, T& val, unsigned cols);
int getnglst(const string& node, vector<string>& nglst);
vector<double> meanprof(const vector<string>& gset, map<string, vector<double> >& exprmap);
template<class T>
double mean(const vector<T>& v);
template<class T>
double stnd(const vector<T>& v, double m);
int get1stcol(const string& f, vector<string>& list);

// Get the node gene list.
int getnglst(const string& f, vector<string>& nglst)
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return 1;
	}

	while(h.good())
	{
		string gene;
		h >> gene;
		h.ignore(1, '\n');
		if(gene == "")
			continue;
		nglst.push_back(gene);
	}
	h.close(); // node gene list read finished!

	return 0;
}

// Get sequences given a set of genes.
int getseq(const string& f, map<string, string>& m)
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return 1;
	}

	while(h.good())
	{
		string gene, seq;
		h >> gene >> seq;
		h.ignore(1, '\n');
		if(gene == "")
			continue;
		str2upper(gene);
		str2upper(seq);
		m[gene] = seq;
	}
	h.close();

	return 0;
}

// Transform a string to upper case.
string& str2upper(string& str)
{
	for(unsigned i = 0; i < str.length(); i++)
		str[i] = toupper(str[i]);

	return str;
}

// Pearson correlation of two vectors x & y.
double Pearcorr(const vector<double>& x, const vector<double>& y)
{
	assert(x.size() == y.size());
	double mx = mean(x);
	double my = mean(y);
	double sx = stnd(x, mx);
	double sy = stnd(y, my);

	double sumcov = 0.0;
	for(unsigned i = 0; i < x.size(); i++)
		sumcov += (x.at(i)-mx)*(y.at(i)-my);

	return sumcov/(x.size()*sx*sy);
}

// Uncentered Pearson correlation.
double Pearcorr0(const vector<double>& x, const vector<double>& y)
{
	assert(x.size() == y.size());
	double ssx = 0.0;
	for(size_t i = 0; i < x.size(); i++)
		ssx += x[i]*x[i];
	double ssy = 0.0;
	for(size_t i = 0; i < y.size(); i++)
		ssy += y[i]*y[i];
		
	double sumcov = 0.0;
	for(size_t i = 0; i < x.size(); i++)
		sumcov += x.at(i)*y.at(i);

	return sumcov/(sqrt(ssx)*sqrt(ssy));
}

// Calculate the mean correlation of a set of genes.
// Tag "uncen" decides whether uncentered Pearson correlation is used.
double meancorr(const vector<string>& gset, map<string, vector<double> >& exprmap, bool uncen)
{
	if(gset.size() == 1)
		return 0.0;

	double mc = 0.0;
	unsigned count = 0;
	for(unsigned i = 0; i < gset.size()-1; i++)
	{
		for(unsigned j = i+1; j < gset.size(); j++)
		{
			if(uncen)
				mc += Pearcorr0(exprmap[gset.at(i)], exprmap[gset.at(j)]);
			else
				mc += Pearcorr(exprmap[gset.at(i)], exprmap[gset.at(j)]);
			count++;
		}
	}

	return mc/count;
}

// Calculate the mean expression profile of a set of genes.
vector<double> meanprof(const vector<string>& gset, map<string, vector<double> >& exprmap)
{
	size_t len = exprmap.begin()->second.size();
	vector<double> mv;
	mv.resize(len, 0.0);
	for(size_t i = 0; i < gset.size(); i++)
	{
		for(size_t j = 0; j < len; j++)
			mv[j] += exprmap[gset[i]][j];
	}
	for(size_t i = 0; i < len; i++)
		mv[i] /= gset.size();

	return mv;
}

// Load all gene expression patterns into memory.
int loadexpr(const string& f, map<string, vector<double> >& exprMap, size_t cols)
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return 1;
	}

	string strLn;
	getline(h, strLn); // get rid of headline;
	while(h.good())
	{
		getline(h, strLn);
		string::size_type tab = strLn.find('\t');
		if(tab != string::npos)
		{
			string id = strLn.substr(0, tab);
			str2upper(id);
			id = id.substr(0, id.find(' '));

			//tab = strLn.find('\t', tab + 1);	// Two head columns: ID & Annotation.
			string sExpr = strLn.substr(tab + 1);
			vector<double> val;
			val.resize(cols);
			str2vec(sExpr, val, (unsigned)cols);
			exprMap[id] = val;
		}
	}
	h.close();
	return 0;
}

// Transform a string to a vector of values.
template<class T>
void str2vec(const string& str, T& val, unsigned cols)
{
	basic_string<char>::size_type tabLoc1 = -1, tabLoc2;
	for(unsigned i = 0; i < cols - 1; i++)
	{
		tabLoc2 = str.find('\t', tabLoc1 + 1);
		if(tabLoc2 - 1 == tabLoc1)
			val[i] = MISSING;
		else
		{
			string vStr = str.substr(tabLoc1 + 1, tabLoc2 - tabLoc1 - 1);
			val[i] = atof(vStr.data());
		}
		tabLoc1 = tabLoc2;
	}
	string vStr = str.substr(tabLoc1 + 1, str.length() - tabLoc1);
	val[cols-1] = atof(vStr.data());
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
	
	return list.size();
}





