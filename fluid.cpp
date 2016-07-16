// fluid.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstring>
#include <sstream>
#include "engine.h"
#include "windows.h"

#define Idx(x, y) (x + N * ( y - 1 ))

#define N 320
#define Tmax 3000

using namespace std;

void initialize(double *pos, double *vel, double T, double X, double Y);
void acceletate(double *acc, double *pos, double X, double Y);
Engine* draw_ini();
void draw(Engine *ep, double *pos, double X, double Y);


int _tmain(int argc, _TCHAR* argv[])
{
	const double X = 20, Y = 20, T = 1;
	const double dt = 0.001;
	const double dt2 = dt * dt / 2;
	const double Ek0 = N * T;

	double pos[N * 2], vel[N * 2], acc[N * 2], newpos[N * 2];
	double Ek = 0;

	initialize(pos, vel, T, X, Y);
	acceletate(acc, pos, X, Y);

	Engine *ep;
	ep = draw_ini();

	for (int t = 0; t < Tmax; t++) {
		for (int i = 0; i < N; i++) {
			newpos[Idx(i, 1)] = pos[Idx(i, 1)] + dt * vel[Idx(i, 1)] + dt2 * acc[Idx(i, 1)];
			newpos[Idx(i, 2)] = pos[Idx(i, 2)] + dt * vel[Idx(i, 2)] + dt2 * acc[Idx(i, 2)];
			pos[Idx(i, 1)] = fmod(newpos[Idx(i, 1)] + X, X);
			pos[Idx(i, 2)] = fmod(newpos[Idx(i, 2)] + Y, Y);
			vel[Idx(i, 1)] += acc[Idx(i, 1)] * dt / 2;
			vel[Idx(i, 2)] += acc[Idx(i, 2)] * dt / 2;
		}
		acceletate(acc, pos, X, Y);
		for (int i = 0; i < N * 2; i++) {
			vel[i] += acc[i] * dt / 2;
			Ek += vel[i] * vel[i] / 2;
		}
		for (int i = 0; i < N * 2; i++)
			vel[i] *= sqrt(Ek0 / Ek);
		Ek = 0;

		draw(ep, pos, X, Y);
		cout << "t = " << t * dt << endl;
	}
		
	getchar();
	engClose(ep);
	return 0;
}

void initialize(double *pos, double *vel, double T, double X, double Y)
{
	double Ek0 = N * T;
	double xa = X / (sqrt(N) + 1), ya = Y / (sqrt(N) + 1);
	double x0 = xa / 2, y0 = ya / 2;
	for (int i = 0; i < N; i++)
	{
		pos[Idx(i,1)] = x0;
		pos[Idx(i,2)] = y0;
		x0 += xa;
		if (x0 > X)  x0 -= X, y0 += ya; 
	}

	srand(time(NULL));
	double Ekk = 0;
	for (int i = 0; i < N * 2; i++)
	{
		vel[i] = rand();
		Ekk += vel[i] * vel[i] / 2;
	}
	for (int i = 0; i < N * 2; i++) vel[i] *= sqrt(Ek0 / Ekk);
}

void acceletate(double *acc, double *pos, double X, double Y)
{
	for (int i = 0; i < N * 2; i++) acc[i] = 0;

	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			double x1 = pos[Idx(i, 1)], y1 = pos[Idx(i, 2)];
			double x2, y2;

			double x20 = pos[Idx(j, 1)];
			double x21 = x20 + X, x22 = x20 - X;
			if (abs(x20 - x1) < __min(abs(x21 - x1), abs(x22 - x1)))
				x2 = x20;
			else if (abs(x21 - x1) < abs(x22 - x1))
				x2 = x21;
			else
				x2 = x22;

			double y20 = pos[Idx(j, 2)];
			double y21 = y20 + Y, y22 = y20 - Y;
			if (abs(y20 - y1) < __min(abs(y21 - y1), abs(y22 - y1)))
				y2 = y20;
			else if (abs(y21 - y1) < abs(y22 - y1))
				y2 = y21;
			else
				y2 = y22;

			double r2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
			double Ft = 12 / pow(r2, 7) - 6 / pow(r2, 4);
			for (int i = 0; i < N; i++) {
				acc[Idx(i, 1)] += Ft * (x1 - x2);
				acc[Idx(i, 2)] += Ft * (y1 - y2);
				acc[Idx(j, 1)] += Ft * (x2 - x1);
				acc[Idx(j, 2)] += Ft * (y2 - y1);
			}
		}
	}
}

Engine* draw_ini()
{
	Engine* ep;
	const int BUFFER_SIZE = 1024;
	char buffer[BUFFER_SIZE];
	if ((ep = engOpen("")) == NULL)
		cout << "Engine Fail" << endl;
	engOutputBuffer(ep, buffer, BUFFER_SIZE);
	cout << "Init Success" << endl;
	return ep;
}

void draw(Engine *ep, double *pos, double X, double Y)
{
	mxArray *x1 = NULL, *y1 = NULL;
	double x[N], y[N];
	for (int i = 0; i < N; i++)
		x[i] = pos[Idx(i, 1)], y[i] = pos[Idx(i, 2)];
	
	x1 = mxCreateDoubleMatrix(N, 1, mxREAL);
	y1 = mxCreateDoubleMatrix(N, 1, mxREAL);

	memcpy((void *)mxGetPr(x1), (void *)x, sizeof(x));
	memcpy((void *)mxGetPr(y1), (void *)y, sizeof(y));

	engPutVariable(ep, "x", x1);
	engPutVariable(ep, "y", y1);

	stringstream sX, sY;
	sX << int(X); sY << int(Y);
	string ssX = sX.str(), ssY = sY.str();
	string s1 = "plot(x, y, 'o'); axis([0 " + ssX + " 0 " + ssY + " ]);";
	const char *plotchar = s1.c_str();
	engEvalString(ep, plotchar);
	engEvalString(ep, "drawnow");
}