#pragma once
#include "Node.h"
#include <math.h>
class Surface
{
public:

	Node *point1;
	Node *point2;

	double ksi[2];
	double eta[2];

	double macierzN[2][4];

	double macierzNN1[4][4]; // jest juz przemnozona razy alfa convection
	double macierzNN2[4][4]; // jest juz przemnozona razy alfa convection
	double macierzSumNN[4][4];

	double det;
	double dlugoscBoku(double x1, double x2, double y1, double y2);
	double edgeLength;

	bool isEdge();

	double N1(double ksi, double eta);
	double N2(double ksi, double eta);
	double N3(double ksi, double eta);
	double N4(double ksi, double eta);


	Surface();
	Surface(Node &n1, Node &n2, double k1, double k2, double e1, double e2);
	~Surface();
	double wekPKrawedz[4];


	Surface & operator=(const Surface &s)
	{
		this->point1 = s.point1;
		this->point2 = s.point2; //tuuuuu reszta przypisania!!!!
		this->ksi[0] = s.ksi[0];
		this->ksi[1] = s.ksi[1];
		this->eta[0] = s.eta[0];
		this->eta[1] = s.eta[1];
		this->edgeLength = s.edgeLength;
		this->det = s.det;
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 4; j++) {
				this->macierzN[i][j] = s.macierzN[i][j];

			}
		}

		
		return *this;
	}
};



inline double Surface::dlugoscBoku(double x1, double x2, double y1, double y2)
{
	if (x1 == x2) 
	{
		return abs(y1 - y2);
	}
	else
	{
		return abs(x1 - x2);
	}

	return 0.0;
}

inline bool Surface::isEdge()
{
	if (this->point1->isOnEdge && this->point2->isOnEdge)
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline Surface::Surface()
{
}

Surface::Surface(Node &n1, Node &n2, double k1, double k2, double e1, double e2)
{
	this->point1 = &n1;
	this->point2 = &n2;
	this->ksi[0] = k1;
	this->ksi[1] = k2;
	this->eta[0] = e1;
	this->eta[1] = e2;
	this->edgeLength = dlugoscBoku(point1->x, point2->x, point1->y, point2->y);
	this->det = edgeLength / 2;
	for (int i = 0; i < 2; i++) {
		this->macierzN[i][0] = N1(ksi[i], eta[i]);
		this->macierzN[i][1] = N2(ksi[i], eta[i]);
		this->macierzN[i][2] = N3(ksi[i], eta[i]);
		this->macierzN[i][3] = N4(ksi[i], eta[i]);
	}
}


Surface::~Surface()
{
}


inline double Surface::N1(double ksi, double eta)
{
	return 0.25 * (1 - ksi) * (1 - eta);
}

inline double Surface::N2(double ksi, double eta)
{
	return 0.25 * (1 + ksi) * (1 - eta);
}

inline double Surface::N3(double ksi, double eta)
{
	return 0.25 * (1 + ksi) * (1 + eta);
}

inline double Surface::N4(double ksi, double eta)
{
	return 0.25 * (1 - ksi) * (1 + eta);
}
