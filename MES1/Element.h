#pragma once

#include "Node.h"
#include "Surface.h"

class Element
{
public:

	int nodesOfElement[4]; // 4-element array of nodes that make an element; inside are numbers of index of nodes

	double ksi[4]; // wspolrzedne calkowania w objetosci
	double eta[4]; // wspolrzedne calkowania w objetosci
	double pochodnePoKsi[4][4];
	double pochodnePoEta[4][4];
	double macierzH[4][4];

	double macierzJ[4][4];

	double macierzDetJ[4];

	double macierzOdwrJ[4][4];

	double pochodnePoX[4][4];

	double pochodnePoY[4][4];

	double macierzTX1[4][4];
	double macierzTX2[4][4];
	double macierzTX3[4][4];
	double macierzTX4[4][4];

	double macierzTY1[4][4];
	double macierzTY2[4][4];
	double macierzTY3[4][4];
	double macierzTY4[4][4];

	double macierzPodcalkowa1[4][4];
	double macierzPodcalkowa2[4][4];
	double macierzPodcalkowa3[4][4];
	double macierzPodcalkowa4[4][4];

	double N1(double ksi, double eta);
	double N2(double ksi, double eta);
	double N3(double ksi, double eta);
	double N4(double ksi, double eta);

	double macierzN[4][4];
	double macierzNN1[4][4];
	double macierzNN2[4][4];
	double macierzNN3[4][4];
	double macierzNN4[4][4];

	double macierzC[4][4];

	double macierzHWarBrzeg[4][4];

	double macierzHOstateczna[4][4];

	Surface *surface = new Surface[4]; // numered oposite to clock moves, from the bottom





	Element();
	~Element();
};



inline double Element::N1(double ksi, double eta)
{
	return 0.25 * (1 - ksi) * (1 - eta);
}

inline double Element::N2(double ksi, double eta)
{
	return 0.25 * (1 + ksi) * (1 - eta);
}

inline double Element::N3(double ksi, double eta)
{
	return 0.25 * (1 + ksi) * (1 + eta);
}

inline double Element::N4(double ksi, double eta)
{
	return 0.25 * (1 - ksi) * (1 + eta);
}

Element::Element()
{
}


Element::~Element()
{
}
