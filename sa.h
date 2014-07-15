/*	sa.h

	Declarations of global variables of Simulated Annealing.
*/

#ifndef SA_H
#define SA_H

// Simulated Annealing.
extern int Repeat;
extern int Iteration;
extern int Changes;
extern double Alpha;
extern double Initemp;
extern double Temp;
extern int chng;	// Counter: changes of BN.
extern int maxpa;
extern int DeterNum;
extern double DeterScor;
extern bool Restag;
extern int Resthrld;
extern int Restarts;
extern int rests;

// Best solution.
struct BSolu;
extern struct BSolu bsolu;

// Tag for using best solution data structure.
extern bool tagbests;

#endif

