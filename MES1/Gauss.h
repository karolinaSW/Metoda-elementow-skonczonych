#pragma once
#include <math.h>
#include <iostream>
using namespace std;

class Gauss
{
public:
	Gauss();
	virtual ~Gauss();
	double * solve(int n, double **arrUnknown, double* arrResult);
};



Gauss::Gauss()
{
}


Gauss::~Gauss()
{
}

inline double * Gauss::solve(int a, double ** arrUnknown, double * vecKnown)
{
	double m, s;
	const int n = a;
	double *vecResult = new double[n];
	double e = pow(10, -12); // 0;

	for (int i = 0; i < n; i++) {
		vecResult[i] = 0;
	}

	double **tabAB = new double*[n];
	for (int i = 0; i < n; i++) {
		tabAB[i] = new double[n + 1];
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tabAB[j][i] = arrUnknown[j][i];
		}
	}

	for (int i = 0; i < n; i++) {
		tabAB[i][n] = vecKnown[i];
	}

	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			if (abs(tabAB[i][i]) < e) {
				cout << "dzielnik rowny 0!!!" << endl;
				break;
			}

			m = -tabAB[j][i] / tabAB[i][i];
			for (int k = 0; k < n + 1; k++)
				tabAB[j][k] += m * tabAB[i][k];
		}
	}

	for (int i = n - 1; i >= 0; i--) {
		s = tabAB[i][n];
		for (int j = n - 1; j >= 0; j--)
			s -= tabAB[i][j] * vecResult[j];

		if (abs(tabAB[i][i]) < e) {
			cout << "dzielnik rowny 0!!!" << endl;
			break;
		}

		vecResult[i] = s / tabAB[i][i];
	}



	return vecResult;
}

