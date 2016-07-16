// quantumwell.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <mkl.h>

using namespace std;

#define NumOfEigs 15

double potential(double x);

int _tmain(int argc, _TCHAR* argv[])
{
	double dx = 0.001, xmin = -10.0, xmax = 10.0;
	const int N = (xmax - xmin) / dx + 1;
	/*
	DSTEMR to find eigenvalues(functions) for symmetic tridiagonal matrix
	============ MAIN OUTPUT =============
	m     # of eigenvalues   INT
	w     eigenvlues         DOUBLE[m]
	z     eigenvectors       DOUBLE[N*m]
	======================================
	*/
	char jobz = 'V';							// 'N' eigval only; 'V' both eigval and eigfun
	char range = 'I';							// 'A' all eigs; 'V' vl < eig <= vu; 'I' #il - #iu eigs 
	
	// define Hamiltonian
	double *d = new double[N];					// diagnol elements
	double *e = new double[N - 1];				// off-diagnol elements
	for (int i = 0; i < N; i++)
		d[i] = 1 / dx / dx + potential(xmin + dx * i);
	for (int i = 0; i < N - 1; i++)
		e[i] = -0.5 / dx / dx;
	// define Hamiltonian

	double vl = -5,								// min eig,    range = 'V'
		vu = 5;								    // max eig,    range = 'V'
	int il = 1,								    // min # eig,  range = 'I'
		iu = NumOfEigs,			    			// max # eig,  range = 'I'
		m,									    // OUTPUT: # of eigs found
	    ldz = N,                                // Leading Dimension of z (length of each eigfun) >=N
		nzc = NumOfEigs;                        // Number of z Counts (more than found eigs) >=m
	double *w = new double[N];                  // OUTPUT: all eigvals
	double *z = new double[N * NumOfEigs];      // OUTPUT: eigfuns, z(1:ldz,1:nzc): [N * m]
	int *isuppz = new int[NumOfEigs * 2];       // nonzero elements in z, >= m*2, not of much use
	int tryrac = 1,								// 0 free accuracy, 1 high accuracy
		lwork = N * 18,                         // size of workspace work
		liwork = N * 10,                        // size of workspace iwork
		info;                                   // 0 SUCCESS, else something wrong
	double *work = new double[lwork];           // workspace? larger than 18*N
	int *iwork = new int[liwork];               // workspace? larger than 10*N
	
	dstemr(&jobz, &range, &N, d, e, &vl, &vu,   // all parameters are pointers
		&il, &iu, &m, w, z, &ldz, &nzc,
		isuppz, &tryrac, work, &lwork,
		iwork, &liwork, &info);

	delete isuppz, work, iwork;
	
	for (int i = 0; i < m; i++)
		cout << w[i] << endl;
	getchar();
}

double potential(double x)
{
	double V;
	V = 2 / (exp(x) + exp(-x));
	V = -V * V;
	return V;
}