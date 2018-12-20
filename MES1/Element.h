#pragma once

#include "Node.h"

class Element
{
public:

	int nodesOfElement[4]; // 4-element array of nodes that make an element; inside are numbers of index of nodes

	double ksi[4]; // wspolrzedne calkowania
	double eta[4]; // wspolrzedne calkowania
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




	Element();
	~Element();
};



Element::Element()
{
}


Element::~Element()
{
}
