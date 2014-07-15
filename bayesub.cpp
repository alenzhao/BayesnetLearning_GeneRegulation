/*	bayesub.cpp

	Definitions of all mathematical and I/O subroutines.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <math.h>
#include "bayesub.h"
#include "globals.h"
#include "sa.h"
#include "fisher2.h"

// Learn Bayesian network - BBNet.
double bbnet(vector<Constraint>& cons, vector<CPTRow>& cpt, const vector<Case>& genlst)
{
	double s = addpres(0, cons, cpt, 1, genlst);	// Add first motif into Bayesian network.
#ifdef VERBOSE
	cout << "Adding motif " << mscor[0].name << endl;
#endif
	for(size_t i = 0; i < mscor.size();)
	{
		// Delete constraint to improve score.
		if(i != 0)	// This step is added to possibly delete constraints before new constraints are added.
			s = delcons(cons, cpt, s, genlst);
		for(set<int>::const_iterator mi = mbnd.begin(); mi != mbnd.end(); mi++)	// Try all rules for current motifs in list.
		{
			Constraint c;
			// Test a new functional depth.
			if(i != 0)
				s = updepth(*mi, cons, cpt, s, genlst);
			// Add constraint: distance to TSS.
			if(rb[0] == '1' && !chkcons(cons, "tss", *mi, -1))
			{
				c.desc = "tss";
				c.motif0 = *mi;
				c.motif1 = -1;
				s = addcons(cons, cpt, c, s, genlst, tss_thrds, ntsst);
			}
			// Add constraint: orientation.
			if(rb[1] == '1' && !chkcons(cons, "orien", *mi, -1))
			{
				c.desc = "orien";
				c.motif0 = *mi;
				c.motif1 = -1;
				int orien[] = {0, 1};
				s = addcons(cons, cpt, c, s, genlst, orien, 2);
			}
			// Add constraint: second copy.
			if(rb[2] == '1' && !chkcons(cons, "sec", *mi, -1))
			{
				c.desc = "sec";
				c.motif0 = *mi;
				c.motif1 = -1;
				int sec[] = {-1};
				s = addcons(cons, cpt, c, s, genlst, sec, 1);
			}
			// Add constraint: distance between any two motifs.
			if(rb[3] == '1')
			{
				for(set<int>::const_iterator mj = mbnd.begin(); mj != mbnd.end(); mj++)
				{
					if(*mi == *mj)
						continue;
					if(!chkcons(cons, "dist", *mi, *mj))
					{
						c.desc = "dist";
						c.motif0 = *mi;
						c.motif1 = *mj;
						s = addcons(cons, cpt, c, s, genlst, dist_thrds, ndistt);
					}
				}
			}
			// Add constraint: order relative to TSS between any two motifs.
			if(rb[4] == '1')
			{
				for(set<int>::const_iterator mj = mbnd.begin(); mj != mbnd.end(); mj++)
				{
					if(*mi == *mj)
						continue;
					if(!chkcons(cons, "order", *mi, *mj))
					{
						c.desc = "order";
						c.motif0 = *mi;
						c.motif1 = *mj;
						int order[] = {0, 1};
						s = addcons(cons, cpt, c, s, genlst, order, 2);
					}
				}
			}
			// Add constraint: looping of any two motifs.
			if(rb[5] == '1')
			{
				for(set<int>::const_iterator mj = mbnd.begin(); mj != mbnd.end(); mj++)
				{
					if(*mi == *mj)
						continue;
					if(!chkcons(cons, "loop", *mi, *mj))
					{
						c.desc = "loop";
						c.motif0 = *mi;
						c.motif1 = *mj;
						s = addcons(cons, cpt, c, s, genlst, loop_thrds, nloopt);
					}
				}
			}
		}
		// Delete constraint to improve score.
		s = delcons(cons, cpt, s, genlst);
		// Add one extra updepth step to find correct depth and distance and to improve score.
		for(set<int>::const_iterator mi = mbnd.begin(); mi != mbnd.end(); mi++)
			s = updepth(*mi, cons, cpt, s, genlst);
		// Add constraint: another new motif. If no improvement, break the loop.
		double s1 = s;
		for(;i < mscor.size()-1 && s1 <= s; i++)	// Keep adding motif until improvement achieved OR the last one.
			s1 = addpres((int)i+1, cons, cpt, s, genlst);
		if(s1 <= s)
			break;
		else
		{
			s = s1;
#ifdef VERBOSE
			cout << "Adding motif " << mscor[(int)i].name << endl;
#endif
		}
	}
	return s;
}


// Learn Bayesian network - GBNet.
double gbnet(vector<Constraint>& cons, vector<CPTRow>& cpt, const vector<Case>& genlst)
{
	double s;	// Bayesian score.
	int rep, iter;	// global iterators for repeat AND iteration.
	s = addpres(0, cons, cpt, 1, genlst, true);	// Add first motif into Bayesian network.
#ifdef VERBOSE
	cout << "Adding motif " << mscor[0].name << endl;
#endif
	for(rep = 0; rep < Repeat; rep++)	// repeat level.
	{
		cout << "**** Running Bayesian network at temperature: " << Temp << " ****" << endl;
		for(iter = 0; iter < Iteration; iter++)	// iteration level.
		{
			if(!chkcons(cons, "pres", 0))
			{
				double s1 = addpres(0, cons, cpt, s, genlst, true);	// Add first motif into Bayesian network.
				if(s1 != s)
				{
#ifdef VERBOSE
					cout << "Adding motif " << mscor[0].name << endl;
#endif
					s = s1;
				}
			}
			for(size_t i = 0; i < mscor.size();)
			{
				// Delete constraint to improve score before considering new constraints.
				// s = delcons(cons, cpt, s, genlst, mbnd);
				for(set<int>::const_iterator mi = mbnd.begin(); mi != mbnd.end(); mi++)	// Try all different functional depths for current motifs.
				{
					Constraint c;
					// Test a new functional depth.
					s = updepth(*mi, cons, cpt, s, genlst, true);
					// Add constraint: distance to TSS.
					if(rb[0] == '1' && !chkcons(cons, "tss", *mi))
					{
						c.desc = "tss";
						c.motif0 = *mi;
						c.motif1 = -1;
						s = addcons(cons, cpt, c, s, genlst, tss_thrds, ntsst, true);
					}
					// Add constraint: orientation.
					if(rb[1] == '1' && !chkcons(cons, "orien", *mi))
					{
						c.desc = "orien";
						c.motif0 = *mi;
						c.motif1 = -1;
						int orien[] = {0, 1};
						s = addcons(cons, cpt, c, s, genlst, orien, 2, true);
					}
					// Add constraint: second copy.
					if(rb[2] == '1' && !chkcons(cons, "sec", *mi))
					{
						c.desc = "sec";
						c.motif0 = *mi;
						c.motif1 = -1;
						int sec[] = {-1};
						s = addcons(cons, cpt, c, s, genlst, sec, 1, true);
					}
					// Add constraint: distance between any two motifs.
					if(rb[3] == '1')
					{
						for(set<int>::const_iterator mj = mbnd.begin(); mj != mbnd.end(); mj++)
						{
							if(*mj == *mi)
								continue;
							if(!chkcons(cons, "dist", *mi, *mj))
							{
								c.desc = "dist";
								c.motif0 = *mi;
								c.motif1 = *mj;
								s = addcons(cons, cpt, c, s, genlst, dist_thrds, ndistt, true);
							}
						}
					}
					// Add constraint: order relative to TSS between any two motifs.
					if(rb[4] == '1')
					{
						for(set<int>::const_iterator mj = mbnd.begin(); mj != mbnd.end(); mj++)
						{
							if(*mj == *mi)
								continue;
							if(!chkcons(cons, "order", *mi, *mj))
							{
								c.desc = "order";
								c.motif0 = *mi;
								c.motif1 = *mj;
								int order[] = {0, 1};
								s = addcons(cons, cpt, c, s, genlst, order, 2, true);
							}
						}
					}
					// Add constraint: looping of any two motifs.
					if(rb[5] == '1')
					{
						for(set<int>::const_iterator mj = mbnd.begin(); mj != mbnd.end(); mj++)
						{
							if(*mj == *mi)
								continue;
							if(!chkcons(cons, "loop", *mi, *mj))
							{
								c.desc = "loop";
								c.motif0 = *mi;
								c.motif1 = *mj;
								s = addcons(cons, cpt, c, s, genlst, loop_thrds, nloopt, true);
							}
						}
					}
				}
				// Delete constraint to improve score.
				s = delcons(cons, cpt, s, genlst);
				// Add constraint: another new motif. If no improvement, break the loop.
				double s1 = s;
				for(;i < mscor.size()-1 && s1 == s; i++)
				{
					if(!chkcons(cons, "pres", (int)i+1))
						s1 = addpres((int)i+1, cons, cpt, s, genlst, true);	// "i" refer to passed motif.
				}
				if(s1 == s)
					break;
				else
				{
					s = s1;
#ifdef VERBOSE
					cout << "Adding motif " << mscor[(int)i].name << endl;	// "i" now refer to the new added motif.
#endif
				}
			}	// Iteration.
			if(Restag)
				s = restart(s, cons, cpt, iter);
			if(chng > Changes || rests > Restarts)	// Required changes have been made. Go to next repeat.
				break;				// OR, restarting reaches maximum number.
		}	// Repeat.
		// Summarize information about this repeat.
		cout << chng << " changes have been made at temperature: " << Temp << endl;
		cout << "After " << iter << " iterations." << endl;
		cout << "And " << rests << " restarts." << endl;

		if(chng < Resthrld && iter == Iteration)	// Temperature is cool now.
			Restag = true;	// Set tag to restart SA if bad condition happens.
		if(chng == 0 && rests == 0)	// No changes(restarts) have been made. Exit.
			break;
		chng = 0;	// Reset counter for changes.
		rests = 0;	// Reset counter for restartings.
		Temp *= Alpha;	// Decrease temperature by rate alpha.
	}	// Whole procedure.

	return s;
}
// *************** Main Routines End Here ******************* //


// Take a bootstrap sample for a vector of objects.
vector<Case> bsamp(const vector<Case>& t, const vector<Case>& b)
{
	vector<Case> s;
	for(size_t i = 0; i < t.size(); i++)
	{
		int index = (int)floor(((double)rand()-1)/RAND_MAX*t.size());
		s.push_back(t[index]);
	}
	for(size_t i = 0; i < b.size(); i++)
	{
		int index = (int)floor(((double)rand()-1)/RAND_MAX*b.size());
		s.push_back(b[index]);
	}

	return s;
}

// Calculate number of cases in each category.
int calcnum(const vector<Case>& v, int c)
{
	int n = 0;
	for(size_t i = 0; i < v.size(); i++)
	{
		if(v[i].label == c)
			n++;
	}
	return n;
}

// Find the index in depth array according to the depth.
int depidx(double depth)
{
	for(int i = 0; i < nfunc; i++)
	{
		if(func_depths[i] - depth < 1e-8 && func_depths[i] - depth > -1e-8)
			return i;
	}
	return -1;
}

// Add a presence node into Bayesian network.
double addpres(int mi, vector<Constraint>& cons, vector<CPTRow>& cpt, double s, const vector<Case>& genlst, bool jump)
{
	Constraint c;
	c.desc = "pres";
	c.motif0 = mi;
	c.motif1 = -1;
	c.para = -1;

	string motif = mscor[mi].name;
#ifdef VERBOSE
	cout << "Considering constraint: pres of " << motif;
#endif
	int didx = -1;	// Depth index of original binding, if exist.
	bool tag = false;	// tag to test whether this motif's binding is in stack.
	if(mbnd.find(mi) == mbnd.end())
	{
		tag = true;
		mbnd.insert(mi);
	}
	else
		didx = depidx(mscor[mi].depth);

	double s0 = 1.0;	// Best Bayesian score.
	vector<Constraint> cons0;	// Best Constraints.
	vector<CPTRow> cpt0;	// Best CPT.
	int didx0;	// Best depth index.
	for(int i = 0; i < nfunc; i++)
	{
		vector<Constraint> cons1 = cons;
		vector<CPTRow> cpt1 = cpt;
		mscor[mi].depth = func_depths[i];
		int pres[] = {-1};
		double s1 = addcons(cons1, cpt1, c, s0, genlst, pres, 1, false);	// Add presence without jumping.
		if(s1 > s0 || s0 == 1)
		{
			s0 = s1;
			cons0 = cons1;
			cpt0 = cpt1;
			didx0 = i;
		}
	}
#ifdef VERBOSE
	cout << " ..." << s0 << "(" << s << ")" << endl;
#endif
	if(s == 1 || s0 > s || (1/Temp*(s0-s) > log10((double)rand()/RAND_MAX) && jump))	// Use temperature to control jumping.
	{
#ifdef VERBOSE
		cout << "Accepting constraint: pres of " << motif << endl;
#endif
		chng++;	// Increase counter if accept presence.
		s = s0;
		cons = cons0;
		cpt = cpt0;
		mscor[mi].depth = func_depths[didx0];
		if(tagbests)
			bestsolu(s, cons, cpt, mbnd, mscor);
	}
	else if(tag)	// Motif's binding wasn't in stack, delete it.
		mbnd.erase(mi);
	else	// Motif's binding was in stack, recover it.
		mscor[mi].depth = func_depths[didx];

	return s;
}

// Update functional depth of one motif to improve score.
double updepth(int mi, vector<Constraint>& cons, vector<CPTRow>& cpt, double s, const vector<Case>& genlst, bool jump)
{
	// Test whether there is one constraint contain motif "mi".
	bool tag = false;
	for(size_t i = 0; i < cons.size(); i++)
	{
		if(mi == cons[i].motif0 || mi == cons[i].motif1)
		{
			tag = true;
			break;
		}
	}
	if(!tag)
		return s;

	string motif = mscor[mi].name;
#ifdef VERBOSE
	cout << "Choosing a new depth for motif " << motif;
#endif
	// Backup the original binding index.
	int didx = depidx(mscor[mi].depth);
	double s0 = 1.0;	// Best Bayesian score.
	vector<Constraint> cons0 = cons;	// Best constraints.
	vector<CPTRow> cpt0;	// Best CPT.
	int didx0 = -1;	// Best depth index.
	// Try all different functional depths, update CPT and score if necessary.
	for(int i = 0; i < nfunc; i++)
	{
		if(i == didx)	// skip original depth.
			continue;
		mscor[mi].depth = func_depths[i];
		for(size_t j = 0; j < cons.size(); j++)	// consider all parameters for constraints that contain motif "mi".
		{
			if(mi != cons[j].motif0 && mi != cons[j].motif1)
				continue;

			vector<int> paraset;	// set up a set of parameters for tuning.
			if(cons[j].desc == "pres" || cons[j].desc == "sec")
				paraset.push_back(-1);
			else if(cons[j].desc == "tss")
			{
				for(int k = 0; k < ntsst; k++)
					paraset.push_back(tss_thrds[k]);
			}
			else if(cons[j].desc == "orien" || cons[j].desc == "order")
			{
				for(int k = 0; k < 2; k++)
					paraset.push_back(k);
			}
			else if(cons[j].desc == "dist")
			{
				for(int k = 0; k < ndistt; k++)
					paraset.push_back(dist_thrds[k]);
			}
			else if(cons[j].desc == "loop")
			{
				for(int k = 0; k < nloopt; k++)
					paraset.push_back(loop_thrds[k]);
			}

			// Try different parameters.
			for(size_t k = 0; k < paraset.size(); k++)
			{
				vector<Constraint> cons1 = cons0;
				cons1[j].para = paraset[k];
				vector<CPTRow> cpt1, ppt1;
				constrcpt(cpt1, ppt1, genlst, cons1);
				double s1;
				if(!itag)
					s1 = score((int)cons1.size(), cpt1, ppt1);
				else
					s1 = iscore((int)cons1.size(), cpt1);
				if(s1 > s0 || s0 == 1)	// Store the best data structures.
				{
					cons0 = cons1;
					cpt0 = cpt1;
					s0 = s1;
					didx0 = i;
				}
			}
		}
	}
#ifdef VERBOSE
	cout << " ..." << s0 << "(" << s << ")" << endl;
#endif	
	if(s0 > s || (1/Temp*(s0-s) > log10((double)rand()/RAND_MAX) && jump))	// Use temperature to control jumping.
	{
#ifdef VERBOSE
		cout << "Accepting depth change: " << func_depths[didx0] << "(" << func_depths[didx] << ")" << endl;
#endif
		s = s0;
		cons = cons0;
		cpt = cpt0;
		mscor[mi].depth = func_depths[didx0];
		chng++;	// Increase counter if accept depth change.
		if(tagbests)
			bestsolu(s, cons, cpt, mbnd, mscor);
	}
	else
		mscor[mi].depth = func_depths[didx];

	return s;
}

// Add one new constraint into Bayesian network and update everything if necessary.
double addcons(vector<Constraint>& cons, vector<CPTRow>& cpt, Constraint c, double s, 
			   const vector<Case>& genlst, const int paraset[], int npara, bool jump)
{
	if(c.desc != "pres")
	{
#ifdef VERBOSE
		cout << "Considering constraint: " << c.desc << " of " << mscor[c.motif0].name;
		if(c.motif1 != -1)
			cout << " and " << mscor[c.motif1].name;
#endif
	}
	if((int)cons.size() >= maxpa)
	{
#ifdef VERBOSE
		cout << " ...Reach maximum number of parents...Skip!" << endl;
#endif
		return s;
	}
	double s0, s1;	// "0" means best; "1" means current.
	vector<Constraint> cons0, cons1;
	vector<CPTRow> cpt0, cpt1;
	s0 = 1.0;
	for(int i = 0; i < npara; i++)
	{
		c.para = paraset[i];
		cons1 = cons;
		cons1.push_back(c);
		vector<CPTRow> ppt1;
		constrcpt(cpt1, ppt1, genlst, cons1);
		if(!itag)
			s1 = score((int)cons1.size(), cpt1, ppt1);
		else
			s1 = iscore((int)cons1.size(), cpt1);
		if(s1 > s0 || s0 == 1)
		{
			s0 = s1;
			cons0 = cons1;
			cpt0 = cpt1;
		}
	}
	if(c.desc != "pres")
	{
#ifdef VERBOSE
		cout << " ..." << s0 << "(" << s << ")" << endl;
#endif
	}
	if(s0 > s || s == 1 || (1/Temp*(s0-s) > log10((double)rand()/RAND_MAX) && jump))	// Use jumping depends on switch.
	{
		s = s0;
		cons = cons0;
		cpt = cpt0;
		if(c.desc != "pres")
		{
#ifdef VERBOSE
			cout << "Accepting constraint: " << c.desc << endl;
#endif
			chng++;	// Increase counter if accept adding constraint.
			if(tagbests)
				bestsolu(s, cons, cpt, mbnd, mscor);
		}
	}
	return s;
}

// Delete constraint to improve score.
double delcons(vector<Constraint>& cons, vector<CPTRow>& cpt, double s, const vector<Case>& genlst)
{
	bool tag = true;
	while(tag)
	{
		tag = false;
		for(size_t i = 0; i < cons.size(); i++)
		{
			if(cons.size() <= 1)	// stop before all constraints are removed.
				break;
#ifdef VERBOSE
			cout << "Deleting constraint " << cons[i].desc << " of " << mscor[cons[i].motif0].name;
			if(cons[i].motif1 != -1)
				cout << " and " << mscor[cons[i].motif1].name;
#endif
			vector<Constraint> cons1 = cons;
			cons1.erase(cons1.begin() + i);
			vector<CPTRow> cpt1, ppt;
			constrcpt(cpt1, ppt, genlst, cons1);
			double s1;
			if(!itag)
				s1 = score((int)cons1.size(), cpt1, ppt);
			else
				s1 = iscore((int)cons1.size(), cpt1);
#ifdef VERBOSE
			cout << " ..." << s1 << "(" << s << ")" << endl;
#endif
			if(s1 > s)	// Deletion is greedy.
			{
#ifdef VERBOSE
				cout << "Accepting deletion" << endl;
#endif
				tag = true;
				cons = cons1;
				cpt = cpt1;
				s = s1;
				chng++;	// Increase counter if accept deletion.
				break;
			}
		}
	}

	// Check whether there is any hanging binding for motifs that no longer exists in parents.
	vector<int> vdel;	// array for motif ids to be deleted.
	for(set<int>::const_iterator mi = mbnd.begin(); mi != mbnd.end(); mi++)
	{
		bool tag = false;	// Assume no longer exist.
		for(size_t j = 0; j < cons.size(); j++)
		{
			if(*mi == cons[j].motif0 || *mi == cons[j].motif1)
			{
				tag = true;
				break;
			}
		}
		if(!tag)	// Delete the motif's binding if it no longer exists in parents.
			vdel.push_back(*mi);
	}
	for(size_t i = 0; i < vdel.size(); i++)
		mbnd.erase(vdel[i]);

	if(tagbests)
		bestsolu(s, cons, cpt, mbnd, mscor);

	return s;
}

// Check whether a constraint has already been added.
bool chkcons(const vector<Constraint>& cons, const string& desc, int motif0, int motif1)
{
	for(size_t i = 0; i < cons.size(); i++)
	{
		if(cons[i].desc == desc && (cons[i].motif0 == motif0 && cons[i].motif1 == motif1 || 
			cons[i].motif0 == motif1 && cons[i].motif1 == motif0))
			return true;
	}
	return false;
}

// Format and output the results of Bayesian network on a cluster.
int outbayes(ofstream& hOut, double s, const vector<Constraint>& cons, const vector<CPTRow>& cpt, const vector<MotifScore>& vms, size_t node, size_t bkg)
{
	// Header information.
	hOut << "************ Bayesian network parameters & results ************" << endl << endl;
	if(tagbests)
	{
		hOut << "Number of repeats: " << Repeat << endl;
		hOut << "Number of iterations: " << Iteration << endl;
		hOut << "Number of required changes: " << Changes << endl;
		hOut << "Temperature changing rate Alpha: " << Alpha << endl;
		hOut << "Initial temperature: " << Initemp << endl << endl;
	}
	hOut << "Candidate motifs: " << motifcand << endl;
	hOut << "Use prior information for preferred motifs? " << prior << endl;
	hOut << "Parameter of prior of network structure: " << logK << endl;
	hOut << endl;
	if(!primo.empty() && prior == 1)
	{
		hOut << "Prior counts for preferred motifs: " << pricnt << endl;
		hOut << "Preferred motifs: " << endl;
		for(set<string>::const_iterator i = primo.begin(); i != primo.end(); i++)
			hOut << *i << endl;
	}
	hOut << endl << "Single node presence score and optimal depth:" << endl;
	outscor(hOut, vms);

	hOut << endl << "Bayesian score of the network: " << s << endl;
	hOut << "Un-penalized Bayesian score: " << s + logK*cons.size() << endl;
	hOut << endl << "Constraints: " << endl;
	for(size_t i = 0; i < cons.size(); i++)
	{
		hOut << i + 1 << ". ";
		outcons(hOut, cons[i], mscor);
	}

	hOut << endl << "Number of genes that satisfy each constraint: " << endl;
	hOut << "\tIn bkg\tIn node\tP-value" << endl;
	vector<CPTRow> nv = ebitcpt(cpt, cons.size());
	for(size_t i = 0; i < nv.size(); i++)
		hOut << i + 1 << "\t" << nv[i].k0 << "\t" << nv[i].k1 << "\t" 
		<< fpval(nv[i].k1, node-nv[i].k1, nv[i].k0, bkg-nv[i].k1) << endl;

	hOut << endl << "Conditional probability table:" << endl;
	hOut << "\tk = 0\tk = 1\tP-value" << endl;
	for(size_t i = 0; i < cpt.size(); i++)
	{
		hOut << fmtbinary((int)i, cons.size()) << "\t" << cpt[i].k0 << "\t" << cpt[i].k1;
		if(i != 0)
			hOut << "\t" << fpval(cpt[i].k1, node-cpt[i].k1, cpt[i].k0, bkg-cpt[i].k0) << endl;
		else
			hOut << endl;
	}

	if(tagbests)
	{
		hOut << endl << "* " << bsolu.s << endl;
		for(size_t i = 0; i < bsolu.cons.size(); i++)
		{
			hOut << "* " << i + 1 << ". ";
			outcons(hOut, bsolu.cons[i], bsolu.mscor);
		}

		hOut << endl << "* \tIn bkg\tIn node" << endl;
		nv = ebitcpt(bsolu.cpt, bsolu.cons.size());
		for(size_t i = 0; i < nv.size(); i++)
			hOut << "* " << i + 1 << "\t" << nv[i].k0 << "\t" << nv[i].k1 << endl;

		hOut << endl << "* \tk = 0\tk = 1" << endl;
		for(size_t i = 0; i < bsolu.cpt.size(); i++)
			hOut << "* " << fmtbinary((int)i, bsolu.cons.size()) << "\t" << bsolu.cpt[i].k0 << "\t" << bsolu.cpt[i].k1 << endl;
	}


	return 0;
}

// The number of genes that satisfy each constraint from CPT.
vector<CPTRow> ebitcpt(const vector<CPTRow>& cpt, size_t nc)
{
	vector<CPTRow> nv;
	for(size_t i = 0; i < nc; i++)
	{
		int mask = (int)pow((double)2, (int)i);
		CPTRow np;
		np.k0 = np.k1 = 0;
		for(size_t j = 0; j < cpt.size(); j++)
		{
			if((mask & j) != 0)
			{
				np.k0 += cpt[j].k0;
				np.k1 += cpt[j].k1;
			}
		}
		nv.push_back(np);
	}
	return nv;
}

// Output one constraint using file handle.
void outcons(ofstream& h, const Constraint& c, const vector<MotifScore>& mscor)
{
	if(c.desc == "pres")
		h << "Presence of " << mscor[c.motif0].name << ":" << mscor[c.motif0].depth << endl;
	else if(c.desc == "tss")
		h << "Distance to TSS of " << mscor[c.motif0].name << ":" << c.para << ", " << mscor[c.motif0].depth << endl;
	else if(c.desc == "orien")
	{
		h << "Orientation of " << mscor[c.motif0].name << ":";
		if(c.para == 0)
			h << "F";
		else if(c.para == 1)
			h << "R";
		h << ", " << mscor[c.motif0].depth << endl;
	}
	else if(c.desc == "sec")
		h << "Second copy of " << mscor[c.motif0].name << ", " << mscor[c.motif0].depth << endl;
	else if(c.desc == "dist")
		h << "Distance between " << mscor[c.motif0].name << " and " << mscor[c.motif1].name << ":" << c.para 
		<< ", (" << mscor[c.motif0].depth << "," << mscor[c.motif1].depth << ")" << endl;
	else if(c.desc == "order")
	{
		if(c.para == 0)
			h << mscor[c.motif0].name << " is before " << mscor[c.motif1].name << ":" << c.para 
			<< ", (" << mscor[c.motif0].depth << "," << mscor[c.motif1].depth << ")" << endl;
		else if(c.para == 1)
			h << mscor[c.motif1].name << " is before " << mscor[c.motif0].name << ":" << c.para 
			<< ", (" << mscor[c.motif1].depth << "," << mscor[c.motif0].depth << ")" << endl;
	}
	else if(c.desc == "loop")
		h << "Looping of " << mscor[c.motif0].name << " and " << mscor[c.motif1].name << ":" << c.para 
		<< ", (" << mscor[c.motif0].depth << "," << mscor[c.motif1].depth << ")" << endl;
}

// Output all motif scores and optimal functional depths to file.
void outscor(ofstream& h, const vector<MotifScore>& mscor)
{
	for(size_t i = 0; i < mscor.size(); i++)
		h << mscor[i].name << "\t" << mscor[i].score << "\t" << mscor[i].depth << endl;
}

// Format a non-negative integer to a string using binary representation.
string fmtbinary(int n, size_t t)
{
	string sb;	// string of bits.
	sb.resize(t);
	for(size_t i = 0; i < t; i++)
	{
		int mask = (int)pow((double)2, (int)i);
		if((n & mask) != 0)
			sb[t-1-i] = '1';
		else
			sb[t-1-i] = '0';
	}
	return sb;
}

// Load motif scores from file and sort them in descending order.
int loadscor(vector<MotifScore>& mscor, const string& s)
{
	ifstream hScor(s.data());
	if(!hScor)
	{
		cerr << "Can't open " << s << endl;
		return 1;
	}

	while(hScor.good())
	{
		string strLn;
		getline(hScor, strLn);
		if(strLn == "")
			continue;
		istringstream strmLn(strLn);
		MotifScore ascor;
		strmLn >> ascor.name >> ascor.score >> ascor.depth;
		if(ascor.name[0] == '*')
			continue;
		mscor.push_back(ascor);
	}
	hScor.close();

	sort(mscor.begin(), mscor.end(), cmp);
	if(mscor.end() - mscor.begin() > motifcand)	// prevent deletion out of range.
		mscor.erase(mscor.begin() + motifcand, mscor.end());

	return 0;
}

// Load gene list from file.
int loadgene(vector<Case>& tlst, vector<Case>& blst, const string& n, const string& b)
{
	// Genes in the cluster: label = 1.
	ifstream hGen(n.data());
	if(!hGen)
	{
		cerr << "Can't open " << n << endl;
		return 1;
	}

	while(hGen.good())
	{
		string strLn;
		getline(hGen, strLn);
		if(strLn == "")
			continue;
		istringstream strmLn(strLn);
		string gene;
		strmLn >> gene;
		str2upper(gene);
		Case c;
		c.name = gene;
		c.label = 1;
		tlst.push_back(c);
	}
	hGen.close();

	// Genes in the background: label = 0.
	ifstream hBkg(b.data());
	if(!hBkg)
	{
		cerr << "Can't open " << b << endl;
		return 1;
	}

	while(hBkg.good())
	{
		string strLn;
		getline(hBkg, strLn);
		if(strLn == "")
			continue;
		istringstream strmLn(strLn);
		string gene;
		strmLn >> gene;
		str2upper(gene);
		Case c;
		c.name = gene;
		c.label = 0;
		blst.push_back(c);
	}
	hBkg.close();

	return 0;
}

// Load motif list from file.
int loadmotif(vector<string>& motiflst, const string& f)
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return 1;
	}
	while(h.good())
	{
		string strLn;
		getline(h, strLn);
		if(strLn == "")
			continue;
		istringstream strmLn(strLn);
		string motif;
		strmLn >> motif;
		motiflst.push_back(motif);
	}
	h.close();

	return 0;
}

// Load one motif's binding.
int loadone(GBMap& onebind, const string& motif, const set<string>& genset, const string& folder)
{
	string fBind = folder + "/" + motif + ".func";
	ifstream hBind(fBind.data());
	if(!hBind)
	{
		cerr << "Can't open " << fBind << endl;
		return 1;
	}
	while(hBind.good())
	{
		string strLn;
		getline(hBind, strLn);
		if(strLn == "")
			continue;
		istringstream gStrm(strLn);
		string gene;
		int nb;
		gStrm >> gene >> nb;
		str2upper(gene);
		if(genset.find(gene) == genset.end())
			continue;
		VGB vgb;
		for(int j = 0; j < nb; j++)
		{
			string site;
			gStrm >> site;
			vgb.e.push_back(extrbnd(site));
		}
		onebind.e[gene] = vgb;
	}
	return 0;
}

// Load all motifs' binding.
int loadbind(MotifMap& allbind, const vector<string>& motiflst, const set<string>& genset, const string& folder)
{
	for(size_t i = 0; i < motiflst.size(); i++)
	{
		const string& motif = motiflst[i];
		//GBMap onemap;	// Gene binding map for one motif only.
		//pair<string, GBMap> elem(motif, onemap);	// An element for motif map.
		//allbind.e.insert(elem);
		if(loadone(allbind.e[motif], motif, genset, folder) != 0)
			return 1;
	}
	return 0;
}

// A version for motif score list.
int loadbind(MotifMap& allbind, const vector<MotifScore>& mscor, const set<string>& genset, const string& folder)
{
	for(size_t i = 0; i < mscor.size(); i++)
	{
		const string& motif = mscor[i].name;
		//GBMap onemap;	// Gene binding map for one motif only.
		//pair<string, GBMap> elem(motif, onemap);	// An element for motif map.
		//allbind.e.insert(elem);
		if(loadone(allbind.e[motif], motif, genset, folder) != 0)
			return 1;
	}
	return 0;
}

// Extract binding of a site from a string.
GBinding extrbnd(const string& s)
{
	GBinding gb;
	string orien = s.substr(0, s.find(','));
	gb.orien = orien[0];
	size_t semi = s.find(',', 2);
	string score = s.substr(2, semi - 2);
	gb.score = atof(score.data());
	string loc = s.substr(semi + 1);
	gb.loc = atoi(loc.data());

	return gb;
}

// Comparing routine for sorting motif scores.
bool cmp(MotifScore s0, MotifScore s1)
{
	// Arrange the TRANSFAC motifs on top if necessary.
	/*if((s0.name[3] != 'X' || !isnum(s0.name[2])) && (s1.name[3] == 'X' && isnum(s1.name[2])))
	return true;
	else if((s0.name[3] == 'X' && isnum(s0.name[2])) && (s1.name[3] != 'X' || !isnum(s1.name[2])))
	return false;
	else*/
	return s0.score > s1.score;
}

// Check whether a character is a decimal number.
bool isnum(char c)
{
	if(c - '0' >=0 && c - '0' <= 9)
		return true;
	else
		return false;
}

// Display the candidate motifs and their scores.
void dispscor(const vector<MotifScore>& mscor)
{
	for(size_t i = 0; i < mscor.size(); i++)
		cout << mscor[i].name << "\t" << mscor[i].score << "\t" << mscor[i].depth << endl;
}

// Test whether a gene satisfies one constraint.
int test(const Constraint& c, const VGB& m0, double d0, const VGB& m1, double d1, int tss)
{
	if(c.desc == "pres")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score >= d0)
				return 1;
		}
	}
	else if(c.desc == "tss")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score >= d0 && abs(m0.e[i].loc-tss) <= c.para)
				return 1;
		}
	}
	else if(c.desc == "orien")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score < d0)
				continue;
			if(m0.e[i].orien == 'F' && c.para == 0)
				return 1;
			else if(m0.e[i].orien == 'R' && c.para == 1)
				return 1;
		}
	}
	else if(c.desc == "sec")
	{
		int count = 0;
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score >= d0)
				count++;
			if(count > 1)
				return 1;
		}
	}
	else if(c.desc == "dist")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score < d0)
				continue;
			for(size_t j = 0; j < m1.e.size(); j++)
			{
				if(m1.e[j].score < d1)
					continue;
				if(abs(m0.e[i].loc - m1.e[j].loc) <= c.para)
					return 1;
			}
		}
	}
	else if(c.desc == "order")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score < d0)
				continue;
			for(size_t j = 0; j < m1.e.size(); j++)
			{
				if(m1.e[j].score < d1)
					continue;
				if(m0.e[i].loc < m1.e[j].loc && c.para == 0)
					return 1;
				else if(m0.e[i].loc > m1.e[j].loc && c.para == 1)
					return 1;
			}
		}
	}
	else if(c.desc == "loop")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score < d0)
				continue;
			for(size_t j = 0; j < m1.e.size(); j++)
			{
				if(m1.e[j].score < d1)
					continue;
				if(abs(m0.e[i].loc - m1.e[j].loc) > c.para)
					return 1;
			}
		}
	}
	else
		return 0;
	return 0;
}

// According to a set of constraints, classify a gene into a category. 
// Different combinations of the constraints are described in the bits of an integer.
int classification(const string& gene, const vector<Constraint>& cons)
{
	int resbits = 0;
	int mask = 1;
	for(size_t i = 0; i < cons.size(); i++)
	{
		int tag;
		const string& mnam0 = mscor[cons[i].motif0].name;
		if(allbind.e[mnam0].e.find(gene) == allbind.e[mnam0].e.end())
			cerr << "Motif " << mnam0 << " binding info absent for gene: " << gene << " assume no binding" << endl;
		const VGB& m0 = allbind.e[mnam0].e[gene];
		double depth0 = mscor[cons[i].motif0].depth;
		if(cons[i].desc == "pres" || cons[i].desc == "tss" || 
			cons[i].desc == "orien" || cons[i].desc == "sec")
		{
			const VGB& m1 = m0;	// m1 is NULL.
			double depth1 = -1;	// depth1 is invalid.
			int tss;
			if(mtss.find(gene) != mtss.end())	// only for tss rule.
				tss = mtss[gene];
			else
				tss = 0;
			tag = test(cons[i], m0, depth0, m1, depth1, tss);
		}
		else if(cons[i].desc == "dist" || cons[i].desc == "order" || 
			cons[i].desc == "loop")
		{
			const string& mnam1 = mscor[cons[i].motif1].name;
			if(allbind.e[mnam1].e.find(gene) == allbind.e[mnam1].e.end())
				cerr << "Motif " << mnam1 << " binding info absent for gene: " << gene << " assume no binding" << endl;
			const VGB& m1 = allbind.e[mnam1].e[gene];
			double depth1 = mscor[cons[i].motif1].depth;
			tag = test(cons[i], m0, depth0, m1, depth1);
		}
		if(tag == 1)	// the gene satisfy the constraint.
			resbits |= mask;	// bit operation to set the correponding bit to 1.
		mask <<= 1;	// shift the mask to the next bit position.
	}

	return resbits;
}

// Construct conditional probability table given gene list, constraints and motif binding.
void constrcpt(vector<CPTRow>& cpt, vector<CPTRow>& ppt, const vector<Case>& genlst, const vector<Constraint>& cons)
{
	if(cons.size() < 1)
	{
		cpt.clear();
		return;
	}
	// Set prior CPT.
	setprior(ppt, cons);
	// Initialize the CPT.
	initcpt(cpt, (size_t)pow((double)2, (int)cons.size()));
	// Classify each gene and increase the corresponding CPT entry by one.
	for(size_t i = 0; i < genlst.size(); i++)
	{
		int tidx = classification(genlst[i].name, cons);
		if(genlst[i].label == 0)
			cpt[tidx].k0++;
		else if(genlst[i].label == 1)
			cpt[tidx].k1++;
	}
}

// Initialize CPT.
void initcpt(vector<CPTRow>& cpt, size_t ns, int val)
{
	cpt.resize(ns);
	for(size_t i = 0; i < cpt.size(); i++)
	{
		cpt[i].k0 = val;
		cpt[i].k1 = val;
	}
}

// Add prior information into CPT.
void setprior(vector<CPTRow>& ppt, const vector<Constraint>& cons)
{
	initcpt(ppt, (size_t)pow((double)2, (int)cons.size()), 1);
	if(prior == 0)
		return;
	int mask = 0;
	for(size_t i = 0; i < cons.size(); i++)	// add prior counts to table entry if corresponds to preferred motif.
	{
		if(cons[i].desc == "pres" && primo.find(mscor[cons[i].motif0].name) != primo.end())
			mask |= (int)pow((double)2, (int)i);
	}
	ppt[mask].k1 += pricnt;
}

// Calculate Bayesian score given CPT and priors.
double score(int np, const vector<CPTRow>& cpt, const vector<CPTRow>& ppt)
{
	if(np < 1 || cpt.empty())
		return 1.0;

	double s = -np*logK;
	for(size_t i = 0; i < cpt.size(); i++)
	{
		s += -logamma(ppt[i].k0 + ppt[i].k1 + cpt[i].k0 + cpt[i].k1) - logamma(ppt[i].k0) - logamma(ppt[i].k1);
		s += logamma(ppt[i].k0 + ppt[i].k1) + logamma(ppt[i].k0 + cpt[i].k0) + logamma(ppt[i].k1 + cpt[i].k1);
	}
	return s;
}

// Calculate Normalized Mutual Information given CPT.
double iscore(int np, const vector<CPTRow>& cpt)
{
	if(np < 1 || cpt.empty())
		return 0.0;

	vector<double> vm;	// vector for motif marginal probability.
	vector<double> ve;	// vector for expression marginal probability.
	vector<vector<double> > ft;	// frequency table.
	vm.resize(cpt.size(), 0);
	ve.resize(2, 0);
	ft.resize(cpt.size());
	for(size_t i = 0; i < cpt.size(); i++)
	{
		vm[i] = cpt[i].k0 + cpt[i].k1 + 2;
		ve[0] += cpt[i].k0 + 1;
		ve[1] += cpt[i].k1 + 1;
		ft[i].push_back(cpt[i].k0 + 1);
		ft[i].push_back(cpt[i].k1 + 1);
	}
	double N = ve[0] + ve[1];	// N is total count = number of genes + priors.

	for(size_t i = 0; i < cpt.size(); i++)	// Normalize to probability.
	{
		vm[i] /= N;
		ft[i][0] /= N;
		ft[i][1] /= N;
	}
	ve[0] /= N;
	ve[1] /= N;

	// Motivated by AIC. Add penalization when the number of free parameters increases.
	double muinfo = -np*logK;	// mutual information between motif and expression.
	for(size_t i = 0; i < vm.size(); i++)
	{
		for(size_t j = 0; j < ve.size(); j++)
			muinfo += ft[i][j]*log2(ft[i][j]/(vm[i]*ve[j]));
	}

	return muinfo;
}


// Transform a string to upper case.
string& str2upper(string& str)
{
	for(size_t i = 0; i < str.length(); i++)
		str[i] = toupper(str[i]);

	return str;
}

// Calculate the log Gamma value given a integer.
double logamma(int x)
{
	if(x <= 2)
		return 0.0;
	double v = 0.0;
	for(int i = 2; i < x; i++)
	{
		v += log10((double)i);
	}
	return v;
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

// Using rules, CPT and binding to predict cases, return prediction results.
Pred predict(const vector<Constraint>& cons, const vector<CPTRow>& cpt, const vector<string>& plst, const vector<string>& nlst)
{
	Pred resu;
	resu.TP = -1;
	resu.FP = -1;
	resu.TN = -1;
	resu.FN = -1;

	int TP = 0;
	for(size_t i = 0; i < plst.size(); i++)
	{
		string gene = plst[i];
		int idx = classification(gene, cons);
		if(cpt[idx].k0 <= cpt[idx].k1)
			TP++;
	}
	resu.TP = TP;
	resu.FN = (int)plst.size() - TP;

	int TN = 0;
	for(size_t i = 0; i < nlst.size(); i++)
	{
		string gene = nlst[i];
		int idx = classification(gene, cons);
		if(cpt[idx].k0 > cpt[idx].k1)
			TN++;
	}
	resu.TN = TN;
	resu.FP = (int)nlst.size() - TN;

	return resu;
}

// Predict each gene's probability of being in this cluster and output the list as in Beer's prediction.
vector<BPred> predict(const vector<Constraint>& cons, const vector<CPTRow>& cpt, const vector<string>& genlst, int label)
{
	vector<BPred> vbpred;
	for(size_t i = 0; i < genlst.size(); i++)
	{
		BPred bpred;
		string gene = genlst[i];
		int idx = classification(gene, cons);
		bpred.prob = (double)cpt[idx].k1/(cpt[idx].k0+cpt[idx].k1);
		bpred.label = label;
		bpred.name = gene;
		vbpred.push_back(bpred);
	}

	return vbpred;
}

// Operator overload for different parameter.
vector<BPred> predict(const vector<Constraint>& cons, const vector<CPTRow>& cpt, const vector<Case>& genlst, int label)
{
	vector<BPred> vbpred;
	for(size_t i = 0; i < genlst.size(); i++)
	{
		BPred bpred;
		string gene = genlst[i].name;
		int idx = classification(gene, cons);
		bpred.prob = (double)cpt[idx].k1/(cpt[idx].k0+cpt[idx].k1);
		bpred.label = label;
		bpred.name = gene;
		vbpred.push_back(bpred);
	}

	return vbpred;
}

// Output prediction results.
int outpred(ofstream& h, Pred d, const string& node, const string& bkg, const string& pos, const string& neg)
{
	h << endl;
	h << "**** Results of prediction ****" << endl;
	h << "Positive training: " << node << endl;
	h << "Negative training: " << bkg << endl;
	h << "Positive testing: " << pos << endl;
	h << "Negative testing: " << neg << endl;
	h << endl;
	h << "\tTrue\tFalse" << endl;
	h << "Positive: " << d.TP << "\t" << d.FP << endl;
	h << "Negative: " << d.TN << "\t" << d.FN << endl;

	return 0;
}

// Output each gene's probability.
int outpred(ofstream& h, const vector<BPred>& bp)
{
	for(size_t i = 0; i < bp.size(); i++)
		h << bp[i].label << '\t' << bp[i].prob << '\t' << bp[i].name << endl;

	return 0;
}

// Restart SA if bad condition happens.
double restart(double s, vector<Constraint>& cons, vector<CPTRow>& cpt, int& iter)
{
	int ng = 0;	// number of genes satisfying constraints.
	for(size_t i = 1; i < cpt.size(); i++)	// ignore the first row which corresponds to no rules.
		ng += cpt[i].k1;
	if(ng < DeterNum || bsolu.s - s > DeterScor)
	{
#ifdef VERBOSE
		cout << "Bad condition happens! Restart SA with best solution..." << endl;
#endif
		s = bsolu.s;
		cons = bsolu.cons;
		cpt = bsolu.cpt;
		mbnd = bsolu.mbnd;
		mscor = bsolu.mscor;
		chng = 0;	// reset BN counter.
		iter = -1;	// reset iteration counter. looper will automatically add one.
		rests++;	// restarting counter add one.
	}
	return s;
}

// Calculate P-value based on Fisher's exact test.
double fpval(double nm, double nn, double bm, double bn)
{
	double contingency_table[4] = {nm, nn, bm, bn};

	// stuff for Fisher's test.
	int nrow = 2;
	int *nrowp = &nrow;
	int ncol = 2;
	int *ncolp = &ncol;
	double expected = -1.0;
	double percnt = 100.0;
	double emin = 0.0;
	double prt = 0.0;
	double *expectedp = &expected;
	double *percntp = &percnt;
	double *eminp = &emin;
	double *prtp = &prt;
	double pvalue = 0.0;
	double *pvaluep = &pvalue;
	int workspace = 300000;
	int *workspacep = &workspace;
	fexact(nrowp, ncolp, contingency_table, nrowp, expectedp, percntp, eminp, prtp, pvaluep, workspacep);

	return pvalue;
}

// Output each gene's TF binding site information.
int outgene(const string& f, const vector<Case>& n, const vector<Case>& b, const vector<Constraint>& cons)
{
	ofstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return 1;
	}

	h << "** Node **" << endl << endl;
	for(size_t i = 0; i < n.size(); i++)
		outbind(h, n[i].name, cons);

	h << endl;
	h << "** Background ** " << endl << endl;
	for(size_t i = 0; i < b.size(); i++)
		outbind(h, b[i].name, cons);

	return 0;
}

// Output one gene's TF binding site information.
void outbind(ofstream& h, const string& gene, const vector<Constraint>& cons)
{
		h << gene << endl;
		for(size_t j = 0; j < cons.size(); j++)
		{
			const string& mnam0 = mscor[cons[j].motif0].name;
			const VGB& m0 = allbind.e[mnam0].e[gene];
			double depth0 = mscor[cons[j].motif0].depth;
			if(cons[j].desc == "pres" || cons[j].desc == "tss" || 
				cons[j].desc == "orien" || cons[j].desc == "sec")
			{
				const VGB& m1 = m0;	// m1 is NULL.
				double depth1 = -1;	// depth1 is invalid.
				h << "constraint " << j+1 << "\t" << binds(cons[j], m0, depth0, m1, depth1) << endl;
			}
			else if(cons[j].desc == "dist" || cons[j].desc == "order")
			{
				const string& mnam1 = mscor[cons[j].motif1].name;
				const VGB& m1 = allbind.e[mnam1].e[gene];
				double depth1 = mscor[cons[j].motif1].depth;
				h << "constraint " << j+1 << "\t" << binds(cons[j], m0, depth0, m1, depth1) << endl;
			}
		}
}

// The binding sites that satisfy one constraint.
string binds(const Constraint& c, const VGB& m0, double d0, const VGB& m1, double d1)
{
	string info = "";
	if(c.desc == "pres")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score >= d0)
				info += fmtsite(m0.e[i]) + "\t";
		}
	}
	else if(c.desc == "tss")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score >= d0 && m0.e[i].loc <= c.para)
				info += fmtsite(m0.e[i]) + "\t";
		}
	}
	else if(c.desc == "orien")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score < d0)
				continue;
			if(m0.e[i].orien == 'F' && c.para == 0)
				info += fmtsite(m0.e[i]) + "\t";
			else if(m0.e[i].orien == 'R' && c.para == 1)
				info += fmtsite(m0.e[i]) + "\t";
		}
	}
	else if(c.desc == "sec")
	{
		int count = 0;
		string info1 = "";
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score >= d0)
			{
				count++;
				info1 += fmtsite(m0.e[i]) + "\t";
			}
		}
		if(count > 1)
			info = info1;
	}
	else if(c.desc == "dist")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score < d0)
				continue;
			for(size_t j = 0; j < m1.e.size(); j++)
			{
				if(m1.e[j].score < d1)
					continue;
				if(abs((m0.e[i].loc - m1.e[j].loc)) <= c.para)
					info += fmtsite(m0.e[i]) + fmtsite(m1.e[j]) + "\t";
			}
		}
	}
	else if(c.desc == "order")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score < d0)
				continue;
			for(size_t j = 0; j < m1.e.size(); j++)
			{
				if(m1.e[j].score < d1)
					continue;
				if(m0.e[i].loc < m1.e[j].loc && c.para == 0)
					info += fmtsite(m0.e[i]) + fmtsite(m1.e[j]) + "\t";
				else if(m0.e[i].loc > m1.e[j].loc && c.para == 1)
					info += fmtsite(m0.e[i]) + fmtsite(m1.e[j]) + "\t";
			}
		}
	}
	else if(c.desc == "loop")
	{
		for(size_t i = 0; i < m0.e.size(); i++)
		{
			if(m0.e[i].score < d0)
				continue;
			for(size_t j = 0; j < m1.e.size(); j++)
			{
				if(m1.e[j].score < d1)
					continue;
				if(abs((m0.e[i].loc - m1.e[j].loc)) > c.para)
					info += fmtsite(m0.e[i]) + fmtsite(m1.e[j]) + "\t";
			}
		}
	}
	else
		return "NONE";
	if(info != "")
		return info;
	else
		return "NONE";
}

// Format a string describing a binding site.
string fmtsite(const GBinding& gb, bool parth)
{
	ostringstream strm;
	if(parth)
		strm << "(" << gb.orien << "," << gb.score << "," << gb.loc << ")";
	else
		strm << gb.orien << "," << gb.score << "," << gb.loc;
	return strm.str();
}

// Load translational/transcriptional start sites.
int loadtss(const string& f, map<string, int>& m)
{
	ifstream h(f.data());
	if(!h)
	{
		cerr << "Can't open " << f << endl;
		return 1;
	}
	while(h.good())
	{
		string strLn;
		getline(h, strLn);
		size_t tab = strLn.find('\t');
		string gene = strLn.substr(0, tab);
		string tss = strLn.substr(tab+1);
		m[gene] = atoi(tss.data());
	}
	
	return 0;
}


// Record the best solution.
inline void bestsolu(double s, const vector<Constraint>& cons, const vector<CPTRow>& cpt, const set<int> mbnd, const vector<MotifScore>& mscor)
{
	if(bsolu.s == 1 || s > bsolu.s)
	{
		bsolu.s = s;
		bsolu.cons = cons;
		bsolu.cpt = cpt;
		bsolu.mbnd = mbnd;
		bsolu.mscor = mscor;
#ifdef VERBOSE
		cout << "Best solution updated!" << endl;
#endif
	}
}

// Logrithm base 2 of value x.
inline double log2(double x)
{
	return log(x)/log((double)2);
}






