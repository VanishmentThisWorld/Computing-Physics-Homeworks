#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

double myfun(double x);
double integrate_rectangle(double x0, double x1, int N);
double integrate_trapezoidal(double x0, double x1, int N);
double integrate_simpson2(double x0, double x1, int N);
double integrate_simpson3(double x0, double x1, int N);
double integrate_boole(double x0, double x1, int N);
double integrate_montecarlo(double x0, double x1, int M);
double integrate_montecarlo_imp(double x0, double x1, int M);
double gaussrand();

#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616

void main()
{
	const int N = 10000000;
	const int M = 1000000;
	const double x0 = 0;
	const double x1 = 100;
	clock_t start_t, finish_t;
	double total_t;

	cout << "**************************************" << endl;
	cout << "* Integration with different methods *" << endl;
	cout << "* Myfunc      y=exp(-x)              *" << endl;
	cout << "* Lower lim   0                      *" << endl;
	cout << "* Upper lim   100                    *" << endl;
	cout << "**************************************" << endl<<endl;

	start_t = clock();
	double S_rec = integrate_rectangle(x0, x1, N);
	finish_t = clock();
	total_t = (finish_t - start_t) * 1000 / CLOCKS_PER_SEC;
	cout.precision(16);
	cout << "Rectangle:   " << S_rec << endl << "Time: " << total_t << endl << endl;

	start_t = clock();
	double S_tra = integrate_trapezoidal(x0, x1, N);
	finish_t = clock();
	total_t = (finish_t - start_t) * 1000 / CLOCKS_PER_SEC;
	cout.precision(16);
	cout << "Trapezoidal: " << S_tra << endl << "Time: " << total_t << endl << endl;

	start_t = clock();
	double S_si2 = integrate_simpson2(x0, x1, N);
	finish_t = clock();
	total_t = (finish_t - start_t) * 1000 / CLOCKS_PER_SEC;
	cout.precision(16);
	cout << "Simpson:     " << S_si2 << endl << "Time: " << total_t << endl << endl;

	start_t = clock();
	double S_si3 = integrate_simpson3(x0, x1, N);
	finish_t = clock();
	total_t = (finish_t - start_t) * 1000 / CLOCKS_PER_SEC;
	cout.precision(16);
	cout << "Simpson 3/8: " << S_si3 << endl << "Time: " << total_t << endl << endl;

	start_t = clock();
	double S_boo = integrate_boole(x0, x1, N);
	finish_t = clock();
	total_t = (finish_t - start_t) * 1000 / CLOCKS_PER_SEC;
	cout.precision(16);
	cout << "Boole:       " << S_boo << endl << "Time: " << total_t << endl << endl;

	start_t = clock();
	double S_mon = integrate_montecarlo(x0, x1, M);
	finish_t = clock();
	total_t = (finish_t - start_t) * 1000 / CLOCKS_PER_SEC;
	cout.precision(16);
	cout << "Monte Carlo: " << S_mon << endl << "Time: " << total_t << endl << endl;

	start_t = clock();
	double S_moi = integrate_montecarlo_imp(x0, x1, M);
	finish_t = clock();
	total_t = (finish_t - start_t) * 1000 / CLOCKS_PER_SEC;
	cout.precision(16);
	cout << "Important:   " << S_moi << endl << "Time: " << total_t << endl << endl;

	getchar();
}

double myfun(double x)
{
	double y = exp(-x);
	return y;
}

double integrate_rectangle(double x0, double x1, int N)
{
	double dx = (x1 - x0) / N;
	double x = x0;
	double S = 0;
	for (int i = 0; i < N; i++)
	{
		S += myfun(x);
		x += dx;
	}
	return S*dx;
}

double integrate_trapezoidal(double x0, double x1, int N)
{
	double dx = (x1 - x0) / N;
	double x = x0;
	double S = 0;
	for (int i = 0; i < N; i++)
	{
		S += myfun(x) + myfun(x+dx);
		x += dx;
	}
	return S*dx/2;
}

double integrate_simpson2(double x0, double x1, int N)
{
	double dx = (x1 - x0) / N;
	double x = x0;
	double S = 0;
	for (int i = 0; i < N / 2; i++)
	{
		S += 2 * myfun(x) + 4 * myfun(x + dx);
		x += 2 * dx;
	}
	S = (S + myfun(x1) - myfun(x0)) * dx / 3;
	return S;
}

double integrate_simpson3(double x0, double x1, int N)
{
	double dx = (x1 - x0) / N;
	double x = x0;
	double S = 0;
	for (int i = 0; i < N / 3; i++)
	{
		S += 2 * myfun(x) + 3 * myfun(x + dx) + 3 * myfun(x + 2 * dx);
		x += 3 * dx;
	}
	S = (S + myfun(x1) - myfun(x0)) * 3 * dx / 8;
	return S;
}

double integrate_boole(double x0, double x1, int N)
{
	double dx = (x1 - x0) / N;
	double x = x0;
	double S = 0;
	for (int i = 0; i < N / 4; i++)
	{
		S += 14 * myfun(x) + 32 * myfun(x + dx) + 12 * myfun(x + 2 * dx) + 32 * myfun(x + 3 * dx);
		x += 4 * dx;
	}
	S = (S + 7 * myfun(x1) - 7 * myfun(x0)) * 2 * dx / 45;
	return S;
}

double integrate_montecarlo(double x0, double x1, int M)
{
	srand(time(NULL));
	double L = (x1 - x0) / RAND_MAX;
	double S = 0;
	for (int i = 0; i < M; i++)
	{
		S += myfun(rand() * L + x0);
	}
	return S * L * RAND_MAX / M;
}

double integrate_montecarlo_imp(double x0, double x1, int M)
{
	srand(time(NULL));
	double sigma = 10;
	double L = (x1 - x0) / RAND_MAX;
	double S = 0;
	double r;
	for (int i = 0; i < M; i++)
	{
		r = abs(gaussrand() * sigma);
		S += myfun(r) * exp(r * r / 2 / sigma / sigma);
	}
	return S * sigma * sqrt(M_PI_2) / M;
}

double gaussrand()
{
	double u1, u2, v1, v2, s, x;
	do
	{
		u1 = (double)rand() / RAND_MAX;
		u2 = (double)rand() / RAND_MAX;
		v1 = 2 * u1 - 1;
		v2 = 2 * u2 - 1;
		s = v1 * v1 + v2 * v2;
	} while (s >= 1 || s == 0);
	x = v1 * sqrt(-2 * log(s) / s);
	return x;
}