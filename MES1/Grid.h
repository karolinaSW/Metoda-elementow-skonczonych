#pragma once

#include "Element.h"
#include "Data.h"
#include "Gauss.h"

using namespace std;

class Grid
{
public:

	int numberOfElements;
	int numberOfNodes;

	Node *arrOfNodes;
	Element *arrOfElements;

	double **gH;
	double **gC;
	double *gP;

	double **gHC;
	double *gPC;


	void generateGrid(Data &d);

	Grid();
	~Grid();
};



void Grid::generateGrid(Data &d)
{
	numberOfElements = (d.nH - 1) * (d.nL - 1);
	numberOfNodes = d.nH * d.nL;

	arrOfNodes = new Node[numberOfNodes];
	arrOfElements = new Element[numberOfElements];

	gH = new double*[numberOfNodes];
	gC = new double*[numberOfNodes];
	gP = new double [numberOfNodes];
	gHC = new double*[numberOfNodes];
	gPC = new double[numberOfNodes];


	for (int i = 0; i < numberOfNodes; i++) {
		gH[i] = new double[numberOfNodes];
		gC[i] = new double[numberOfNodes];
		gHC[i] = new double[numberOfNodes];
	}

	for (int i = 0; i < numberOfNodes; i++) {
		for (int j = 0; j < numberOfNodes; j++) {
			gH[i][j] = 0;
			gC[i][j] = 0;
			gHC[i][j] = 0;
		}
		gP[i] = 0.0;
	}


	//----------- ta pêtla jedzie po wszystkich elementach siatki
	int k = 0;
	for (int j = 0; j < (d.nL - 1); j++)
	{
		for (int i = 0; i < (d.nH - 1); i++)
		{

			this->arrOfElements[k].nodesOfElement[0] = k + j;
			this->arrOfElements[k].nodesOfElement[1] = k + d.nH + j;
			this->arrOfElements[k].nodesOfElement[2] = k + d.nH + 1 + j;
			this->arrOfElements[k].nodesOfElement[3] = k + 1 + j;
			

			k++;
		}

	}
	//----------- </> ta pêtla jedzie po wszystkich elementach siatki

	/*
	int a, b, c, e;
	a = this->arrOfElements[0].nodesOfElement[0];
	//b = this->arrOfElements[(d.nL - 2) * (d.nH - 1)].nodesOfElement[1];
	b = this->arrOfElements[((d.nH - 1) * (d.nL - 1)) - 1 - (d.nH - 2)].nodesOfElement[1]; // bo 28 (w numeracji od 1 a mamy numeracjê od zera w tablicy)
	//c = this->arrOfElements[(d.nL - 1) * (d.nH - 1) - 1].nodesOfElement[2];
	c = this->arrOfElements[((d.nH - 1) * (d.nL - 1)) - 1].nodesOfElement[2];
	e = this->arrOfElements[d.nH - 2].nodesOfElement[3];
	cout << "Test of nodes's index: " << a << " - " << b << " - " << c << " - " << e << endl << endl;
	cout << " H= " << d.H << " L = " << d.L << endl;
	
	
	a = this->arrOfElements[0].nodesOfElement[0];
	b = this->arrOfElements[(d.nL - 2) * (d.nH - 1)].nodesOfElement[1];
	c = this->arrOfElements[(d.nL - 1) * (d.nH - 1) - 1].nodesOfElement[2];
	e = this->arrOfElements[d.nH - 2].nodesOfElement[3];
	cout << "Test of nodes's index: " << a << " - " << b << " - " << c << " - " << e << endl << endl;
	cout << " H= " << d.H << " L = " << d.L << endl;
	*/

	double deltaX = (d.L / (d.nL - 1));
	double deltaY = (d.H / (d.nH - 1));
	//cout << " deltax= " << deltaX << " deltaY = " << deltaY << endl;
	k = 0;
	for (int j = 0; j < d.nL; j++)
	{
		for (int i = 0; i < d.nH; i++)
		{
			if ( (j == 0) || (j == (d.nL - 1)) || (i == 0) || (i == (d.nH - 1)) )
			{
				this->arrOfNodes[k].isOnEdge = true;
			}

			this->arrOfNodes[k].x = j * deltaX;
			this->arrOfNodes[k].y = i * deltaY;


			k++;

		}

	}
	/*double x = this->arrOfNodes[d.nH * (d.nL - 1)].x;
	double y = this->arrOfNodes[numberOfNodes - 1].y;
	cout << "Test of nodes's coordinates: \t x = " << x << "\t y = " << y << endl << endl;
	

	cout << endl << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements[0].nodesOfElement[i] << "\t";
	}
	cout << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements[12].nodesOfElement[i] << "\t";
	}
	cout << endl;

	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements[17].nodesOfElement[i] << "\t";
	}
	cout << endl;

	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements[5].nodesOfElement[i] << "\t";
	}
	cout << endl;
	*/

	//cout << endl << endl << endl << " ################################################################" << endl;
	

	//----tutaj test dla danych z zajec:
	/*
	this->arrOfNodes[this->arrOfElements[6].nodesOfElement[0]].x = 0.0; //po konkretnym elemencie wyszukujê numery wêz³ów w tym elemencie i dajê je jako index do tablicy ogólmej wêz³ów
	this->arrOfNodes[this->arrOfElements[6].nodesOfElement[1]].x = 0.025;
	this->arrOfNodes[this->arrOfElements[6].nodesOfElement[2]].x = 0.025;
	this->arrOfNodes[this->arrOfElements[6].nodesOfElement[3]].x = 0.0;

	this->arrOfNodes[this->arrOfElements[6].nodesOfElement[0]].y = 0.0;
	this->arrOfNodes[this->arrOfElements[6].nodesOfElement[1]].y = 0.0;
	this->arrOfNodes[this->arrOfElements[6].nodesOfElement[2]].y = 0.025;
	this->arrOfNodes[this->arrOfElements[6].nodesOfElement[3]].y = 0.025;
	*/
	
	
	//----------- ta pêtla jedzie po wszystkich elementach siatki
	k = 0;
	for (int l = 0; l < (d.nL - 1); l++)
	{
		for (int m = 0; m < (d.nH - 1); m++)
		{
			cout << "element " << k << endl << endl;

			this->arrOfElements[k].eta[0] = -(1 / sqrt(3));
			this->arrOfElements[k].eta[1] = -(1 / sqrt(3));
			this->arrOfElements[k].eta[2] =  (1 / sqrt(3));
			this->arrOfElements[k].eta[3] =  (1 / sqrt(3));

			this->arrOfElements[k].ksi[0] = -(1 / sqrt(3));
			this->arrOfElements[k].ksi[1] =  (1 / sqrt(3));
			this->arrOfElements[k].ksi[2] =  (1 / sqrt(3));
			this->arrOfElements[k].ksi[3] = -(1 / sqrt(3));

			//macierz dN1/dKsi  i  dN1/dEta   dla czterech punktow calkowania
			for (int i = 0; i < 4; i++) // i determinuje punkt calkowania
			{
				//this->arrOfElements->pochodnePoKsi[i][0] = (-0.25 * (1 - this->arrOfElements[k].eta[i]));
				//TU KOMBINOWANE
				this->arrOfElements[k].pochodnePoKsi[i][0] = (-0.25 * (1 - this->arrOfElements[k].eta[i])); // 0,1,2,3 - nr funkcji ksztaltu
				this->arrOfElements[k].pochodnePoKsi[i][1] = ( 0.25 * (1 - this->arrOfElements[k].eta[i]));
				this->arrOfElements[k].pochodnePoKsi[i][2] = ( 0.25 * (1 + this->arrOfElements[k].eta[i]));
				this->arrOfElements[k].pochodnePoKsi[i][3] = (-0.25 * (1 + this->arrOfElements[k].eta[i]));

				this->arrOfElements[k].pochodnePoEta[i][0] = (-0.25 * (1 - this->arrOfElements[k].ksi[i]));
				this->arrOfElements[k].pochodnePoEta[i][1] = (-0.25 * (1 + this->arrOfElements[k].ksi[i]));
				this->arrOfElements[k].pochodnePoEta[i][2] = ( 0.25 * (1 + this->arrOfElements[k].ksi[i]));
				this->arrOfElements[k].pochodnePoEta[i][3] = ( 0.25 * (1 - this->arrOfElements[k].ksi[i]));
			}
			///----------------------------------------TUUUUUUU skoñczy³am ----------------------------------------------
			///-----------------------------------------------------------------
			///--------------------------------------------------------------------
			//macierze J dla 4 punktow calkowania
			for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2; 0123 - pe³ny jakobian dla jednego punktu
			{
				//this->arrOfElements[k].macierzJ[i][0] = ((this->arrOfElements[k].pochodnePoKsi[i][0] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[0]].x) + (this->arrOfElements->pochodnePoKsi[i][1] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[1]].x) + (this->arrOfElements->pochodnePoKsi[i][2] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[2]].x) + (this->arrOfElements->pochodnePoKsi[i][3] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[3]].x));
				this->arrOfElements[k].macierzJ[i][0] = ((this->arrOfElements[k].pochodnePoKsi[i][0] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[0]].x) + (this->arrOfElements[k].pochodnePoKsi[i][1] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[1]].x) + (this->arrOfElements[k].pochodnePoKsi[i][2] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[2]].x) + (this->arrOfElements[k].pochodnePoKsi[i][3] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[3]].x));
				this->arrOfElements[k].macierzJ[i][1] = ((this->arrOfElements[k].pochodnePoKsi[i][0] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[0]].y) + (this->arrOfElements[k].pochodnePoKsi[i][1] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[1]].y) + (this->arrOfElements[k].pochodnePoKsi[i][2] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[2]].y) + (this->arrOfElements[k].pochodnePoKsi[i][3] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[3]].y));
				this->arrOfElements[k].macierzJ[i][2] = ((this->arrOfElements[k].pochodnePoEta[i][0] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[0]].x) + (this->arrOfElements[k].pochodnePoEta[i][1] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[1]].x) + (this->arrOfElements[k].pochodnePoEta[i][2] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[2]].x) + (this->arrOfElements[k].pochodnePoEta[i][3] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[3]].x));
				this->arrOfElements[k].macierzJ[i][3] = ((this->arrOfElements[k].pochodnePoEta[i][0] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[0]].y) + (this->arrOfElements[k].pochodnePoEta[i][1] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[1]].y) + (this->arrOfElements[k].pochodnePoEta[i][2] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[2]].y) + (this->arrOfElements[k].pochodnePoEta[i][3] * this->arrOfNodes[this->arrOfElements[k].nodesOfElement[3]].y));
			}

			//macierze det_J dla 4 punktow calkowania
			for (int i = 0; i < 4; i++)
			{
				this->arrOfElements[k].macierzDetJ[i] = (this->arrOfElements[k].macierzJ[i][0] * this->arrOfElements[k].macierzJ[i][3]) - (this->arrOfElements[k].macierzJ[i][1] * this->arrOfElements[k].macierzJ[i][2]);
			}

			//macierze odwrotna do J dla 4 punktow calkowania
			for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
			{
				this->arrOfElements[k].macierzOdwrJ[i][0] =  (this->arrOfElements[k].macierzJ[i][3] / this->arrOfElements[k].macierzDetJ[i]);
				this->arrOfElements[k].macierzOdwrJ[i][1] = -(this->arrOfElements[k].macierzJ[i][1] / this->arrOfElements[k].macierzDetJ[i]);
				this->arrOfElements[k].macierzOdwrJ[i][2] = -(this->arrOfElements[k].macierzJ[i][2] / this->arrOfElements[k].macierzDetJ[i]);
				this->arrOfElements[k].macierzOdwrJ[i][3] =  (this->arrOfElements[k].macierzJ[i][0] / this->arrOfElements[k].macierzDetJ[i]);
			}

			//macierz pochodnych po x i po y
			for (int i = 0; i < 4; i++)
			{
				this->arrOfElements[k].pochodnePoX[0][i] = (this->arrOfElements[k].macierzOdwrJ[i][0] * this->arrOfElements[k].pochodnePoKsi[i][0]) + (this->arrOfElements[k].macierzOdwrJ[i][1] * this->arrOfElements[k].pochodnePoEta[i][0]);
				this->arrOfElements[k].pochodnePoX[1][i] = (this->arrOfElements[k].macierzOdwrJ[i][0] * this->arrOfElements[k].pochodnePoKsi[i][1]) + (this->arrOfElements[k].macierzOdwrJ[i][1] * this->arrOfElements[k].pochodnePoEta[i][1]);
				this->arrOfElements[k].pochodnePoX[2][i] = (this->arrOfElements[k].macierzOdwrJ[i][0] * this->arrOfElements[k].pochodnePoKsi[i][2]) + (this->arrOfElements[k].macierzOdwrJ[i][1] * this->arrOfElements[k].pochodnePoEta[i][2]);
				this->arrOfElements[k].pochodnePoX[3][i] = (this->arrOfElements[k].macierzOdwrJ[i][0] * this->arrOfElements[k].pochodnePoKsi[i][3]) + (this->arrOfElements[k].macierzOdwrJ[i][1] * this->arrOfElements[k].pochodnePoEta[i][3]);

				this->arrOfElements[k].pochodnePoY[0][i] = (this->arrOfElements[k].macierzOdwrJ[i][2] * this->arrOfElements[k].pochodnePoKsi[i][0]) + (this->arrOfElements[k].macierzOdwrJ[i][3] * this->arrOfElements[k].pochodnePoEta[i][0]);
				this->arrOfElements[k].pochodnePoY[1][i] = (this->arrOfElements[k].macierzOdwrJ[i][2] * this->arrOfElements[k].pochodnePoKsi[i][1]) + (this->arrOfElements[k].macierzOdwrJ[i][3] * this->arrOfElements[k].pochodnePoEta[i][1]);
				this->arrOfElements[k].pochodnePoY[2][i] = (this->arrOfElements[k].macierzOdwrJ[i][2] * this->arrOfElements[k].pochodnePoKsi[i][2]) + (this->arrOfElements[k].macierzOdwrJ[i][3] * this->arrOfElements[k].pochodnePoEta[i][2]);
				this->arrOfElements[k].pochodnePoY[3][i] = (this->arrOfElements[k].macierzOdwrJ[i][2] * this->arrOfElements[k].pochodnePoKsi[i][3]) + (this->arrOfElements[k].macierzOdwrJ[i][3] * this->arrOfElements[k].pochodnePoEta[i][3]);
			}

			//###
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzTX1[i][j] = this->arrOfElements[k].macierzDetJ[0] * this->arrOfElements[k].pochodnePoX[i][0] * this->arrOfElements[k].pochodnePoX[j][0];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzTX2[i][j] = this->arrOfElements[k].macierzDetJ[1] * this->arrOfElements[k].pochodnePoX[i][1] * this->arrOfElements[k].pochodnePoX[j][1];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzTX3[i][j] = this->arrOfElements[k].macierzDetJ[2] * this->arrOfElements[k].pochodnePoX[i][2] * this->arrOfElements[k].pochodnePoX[j][2];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzTX4[i][j] = this->arrOfElements[k].macierzDetJ[3] * this->arrOfElements[k].pochodnePoX[i][3] * this->arrOfElements[k].pochodnePoX[j][3];
				}
			}


			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzTY1[i][j] = this->arrOfElements[k].macierzDetJ[0] * this->arrOfElements[k].pochodnePoY[i][0] * this->arrOfElements[k].pochodnePoY[j][0];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzTY2[i][j] = this->arrOfElements[k].macierzDetJ[1] * this->arrOfElements[k].pochodnePoY[i][1] * this->arrOfElements[k].pochodnePoY[j][1];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzTY3[i][j] = this->arrOfElements[k].macierzDetJ[2] * this->arrOfElements[k].pochodnePoY[i][2] * this->arrOfElements[k].pochodnePoY[j][2];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzTY4[i][j] = this->arrOfElements[k].macierzDetJ[3] * this->arrOfElements[k].pochodnePoY[i][3] * this->arrOfElements[k].pochodnePoY[j][3];
				}
			}

			//##

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzPodcalkowa1[i][j] = d.K * (this->arrOfElements[k].macierzTX1[i][j] + this->arrOfElements[k].macierzTY1[i][j]);
					this->arrOfElements[k].macierzPodcalkowa2[i][j] = d.K * (this->arrOfElements[k].macierzTX2[i][j] + this->arrOfElements[k].macierzTY2[i][j]);
					this->arrOfElements[k].macierzPodcalkowa3[i][j] = d.K * (this->arrOfElements[k].macierzTX3[i][j] + this->arrOfElements[k].macierzTY3[i][j]);
					this->arrOfElements[k].macierzPodcalkowa4[i][j] = d.K * (this->arrOfElements[k].macierzTX4[i][j] + this->arrOfElements[k].macierzTY4[i][j]);

				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzH[i][j] = (this->arrOfElements[k].macierzPodcalkowa1[i][j] + this->arrOfElements[k].macierzPodcalkowa2[i][j] + this->arrOfElements[k].macierzPodcalkowa3[i][j] + this->arrOfElements[k].macierzPodcalkowa4[i][j]);
				}
			}

			// ----------------------- Macierz C

			for (int i = 0; i < 4; i++) // i - punkt calkowania, 0123 - funkcja ksztaltu
			{
				this->arrOfElements[k].macierzN[i][0] = this->arrOfElements[k].N1(this->arrOfElements[k].ksi[i], this->arrOfElements[k].eta[i]);
				this->arrOfElements[k].macierzN[i][1] = this->arrOfElements[k].N2(this->arrOfElements[k].ksi[i], this->arrOfElements[k].eta[i]);
				this->arrOfElements[k].macierzN[i][2] = this->arrOfElements[k].N3(this->arrOfElements[k].ksi[i], this->arrOfElements[k].eta[i]);
				this->arrOfElements[k].macierzN[i][3] = this->arrOfElements[k].N4(this->arrOfElements[k].ksi[i], this->arrOfElements[k].eta[i]);
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					//this->arrOfElements[k].macierzNN1[i][j] = (this->arrOfElements[k].N1(this->arrOfElements[k].ksi[i], this->arrOfElements[k].eta[i]) * this->arrOfElements[k].N1(this->arrOfElements[k].ksi[j], this->arrOfElements[k].eta[j]));
					this->arrOfElements[k].macierzNN1[i][j] = this->arrOfElements[k].macierzN[0][i] * this->arrOfElements[k].macierzN[0][j] * this->arrOfElements[k].macierzDetJ[0] * d.c * d.ro;
					this->arrOfElements[k].macierzNN2[i][j] = this->arrOfElements[k].macierzN[1][i] * this->arrOfElements[k].macierzN[1][j] * this->arrOfElements[k].macierzDetJ[1] * d.c * d.ro;
					this->arrOfElements[k].macierzNN3[i][j] = this->arrOfElements[k].macierzN[2][i] * this->arrOfElements[k].macierzN[2][j] * this->arrOfElements[k].macierzDetJ[2] * d.c * d.ro;
					this->arrOfElements[k].macierzNN4[i][j] = this->arrOfElements[k].macierzN[3][i] * this->arrOfElements[k].macierzN[3][j] * this->arrOfElements[k].macierzDetJ[3] * d.c * d.ro;
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzC[i][j] = this->arrOfElements[k].macierzNN1[i][j] + this->arrOfElements[k].macierzNN2[i][j] + this->arrOfElements[k].macierzNN3[i][j] + this->arrOfElements[k].macierzNN4[i][j];
				}
			}
			//----------------------</> Macierz C

			//-----------------------------------------------------plaszczyzny


			// nadanie plaszczyznom wezlow i wartosci ksi eta
			Surface a(this->arrOfNodes[k + l], this->arrOfNodes[k + d.nH + l], -(1/sqrt(3)), (1 / sqrt(3)), -1, -1);
			this->arrOfElements[k].surface[0] = a;
			Surface b(this->arrOfNodes[k + d.nH + l], this->arrOfNodes[k + d.nH + 1 + l], 1, 1, -(1 / sqrt(3)), (1 / sqrt(3)));
			this->arrOfElements[k].surface[1] = b;
			Surface c(this->arrOfNodes[k + d.nH + 1 + l], this->arrOfNodes[k + 1 + l], (1 / sqrt(3)), -(1 / sqrt(3)), 1, 1);
			this->arrOfElements[k].surface[2] = c;
			Surface e(this->arrOfNodes[k + 1 + l], this->arrOfNodes[k + l], -1, -1, (1 / sqrt(3)), -(1 / sqrt(3)));
			this->arrOfElements[k].surface[3] = e;


			//obliczanie macierzy z wartosciami funkcji ksztaltu dla kazdego punktu 
			for (int j = 0; j < 4; j++)
			{
				for (int i = 0; i < 2; i++) // 0123 - nr funkcji ksztaltu, i - punkt calkowania
				{
					this->arrOfElements[k].surface[j].macierzN[i][0] = this->arrOfElements[k].N1(this->arrOfElements[k].surface[j].ksi[i], this->arrOfElements[k].surface[j].eta[i]);
					this->arrOfElements[k].surface[j].macierzN[i][1] = this->arrOfElements[k].N2(this->arrOfElements[k].surface[j].ksi[i], this->arrOfElements[k].surface[j].eta[i]);
					this->arrOfElements[k].surface[j].macierzN[i][2] = this->arrOfElements[k].N3(this->arrOfElements[k].surface[j].ksi[i], this->arrOfElements[k].surface[j].eta[i]);
					this->arrOfElements[k].surface[j].macierzN[i][3] = this->arrOfElements[k].N4(this->arrOfElements[k].surface[j].ksi[i], this->arrOfElements[k].surface[j].eta[i]);
				}
			}

			for (int j = 0; j < 4; j++)
			{
				for (int i = 0; i < 4; i++)
				{
					for (int h = 0; h < 4; h++) 
					{
						this->arrOfElements[k].surface[j].macierzNN1[i][h] = this->arrOfElements[k].surface[j].macierzN[0][i] * this->arrOfElements[k].surface[j].macierzN[0][h] * d.alfa;
						this->arrOfElements[k].surface[j].macierzNN2[i][h] = this->arrOfElements[k].surface[j].macierzN[1][i] * this->arrOfElements[k].surface[j].macierzN[1][h] * d.alfa;

					}
				}
			}

			for (int j = 0; j < 4; j++)
			{
				for (int i = 0; i < 4; i++)
				{
					for (int h = 0; h < 4; h++)
					{
						this->arrOfElements[k].surface[j].macierzSumNN[i][h] = (this->arrOfElements[k].surface[j].macierzNN1[i][h] + this->arrOfElements[k].surface[j].macierzNN2[i][h]) * this->arrOfElements[k].surface[j].det;
					}
				}
			}

			///--------------------ajajajaja for test only
			//cout << "elem " << k << endl;

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{ 
					this->arrOfElements[k].macierzHWarBrzeg[i][j] = (this->arrOfElements[k].surface[0].isEdge() * this->arrOfElements[k].surface[0].macierzSumNN[i][j]) + (this->arrOfElements[k].surface[1].isEdge()*this->arrOfElements[k].surface[1].macierzSumNN[i][j]) + (this->arrOfElements[k].surface[2].isEdge()*this->arrOfElements[k].surface[2].macierzSumNN[i][j]) + (this->arrOfElements[k].surface[3].isEdge() * this->arrOfElements[k].surface[3].macierzSumNN[i][j]);
				}
			}

			// ---------------------^^^^^^^^^^^^^^^^^^^^^^^^------------------ macierz H ostateczna

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements[k].macierzHOstateczna[i][j] = this->arrOfElements[k].macierzH[i][j] +this->arrOfElements[k].macierzHWarBrzeg[i][j];

					//cout << this->arrOfElements[k].macierzHWarBrzeg[i][j] << "\t";
				}
				//cout << endl;
			}
			/*
			cout << "N: 00 01.../ 10 11..." << endl;
			for (int i = 0; i < 4; i++) {
				cout << this->arrOfElements[k].surface[i].macierzN[0][0] << "  " << this->arrOfElements[k].surface[i].macierzN[0][1] << "   " << this->arrOfElements[k].surface[i].macierzN[0][2] << "   " << this->arrOfElements[k].surface[i].macierzN[0][3] << endl;
				cout << this->arrOfElements[k].surface[i].macierzN[1][0] << "  " << this->arrOfElements[k].surface[i].macierzN[1][1] << "   " << this->arrOfElements[k].surface[i].macierzN[1][2] << "   " << this->arrOfElements[k].surface[i].macierzN[1][3] << endl;
				cout << endl;
			}*/

			for (int i = 0; i < 4; i++) {
				//cout << this->arrOfElements[k].surface[i].isEdge() << "    ";
			}
			cout << endl;


			//cout << "peeeeeeeeee kr" << endl << endl;
			//--------------------------------------------------</>wektor P
			for (int i = 0; i < 4; i++) 
			{
				this->arrOfElements[k].surface[i].wekPKrawedz[0] = - (((this->arrOfElements[k].surface[i].macierzN[0][0] * d.alfa * d.tempAmbient) + (this->arrOfElements[k].surface[i].macierzN[1][0] * d.alfa * d.tempAmbient)) * this->arrOfElements[k].surface[i].det);
				this->arrOfElements[k].surface[i].wekPKrawedz[1] = - (((this->arrOfElements[k].surface[i].macierzN[0][1] * d.alfa * d.tempAmbient) + (this->arrOfElements[k].surface[i].macierzN[1][1] * d.alfa * d.tempAmbient)) * this->arrOfElements[k].surface[i].det);
				this->arrOfElements[k].surface[i].wekPKrawedz[2] = - (((this->arrOfElements[k].surface[i].macierzN[0][2] * d.alfa * d.tempAmbient) + (this->arrOfElements[k].surface[i].macierzN[1][2] * d.alfa * d.tempAmbient)) * this->arrOfElements[k].surface[i].det);
				this->arrOfElements[k].surface[i].wekPKrawedz[3] = - (((this->arrOfElements[k].surface[i].macierzN[0][3] * d.alfa * d.tempAmbient) + (this->arrOfElements[k].surface[i].macierzN[1][3] * d.alfa * d.tempAmbient)) * this->arrOfElements[k].surface[i].det);

				//cout << this->arrOfElements[k].surface[i].wekPKrawedz[0] << "\t" << this->arrOfElements[k].surface[i].wekPKrawedz[1] << "\t" << this->arrOfElements[k].surface[i].wekPKrawedz[2] << "\t" << this->arrOfElements[k].surface[i].wekPKrawedz[3] << "\t" << endl;

			}
			//cout << endl << endl << "caly wek P" << endl;
			for (int i = 0; i < 4; i++) {
				this->arrOfElements[k].wektorP[i] = ((this->arrOfElements[k].surface[0].wekPKrawedz[i] * this->arrOfElements[k].surface[0].isEdge())+(this->arrOfElements[k].surface[1].wekPKrawedz[i] * this->arrOfElements[k].surface[1].isEdge())+(this->arrOfElements[k].surface[2].wekPKrawedz[i] * this->arrOfElements[k].surface[2].isEdge())+(this->arrOfElements[k].surface[3].wekPKrawedz[i] * this->arrOfElements[k].surface[3].isEdge()));
				//cout << this->arrOfElements[k].wektorP[i] << "    ";
			}
			//cout << endl << endl;

			//------------------------------------------------</> wektor P

			
			k++;// dont touch it!!! iterator elementu
		}
	}
	//----------- </> ta pêtla jedzie po wszystkich elementach siatki











	//---------------------------------------------------------------------------------------------------------------------------------------
	//-----------------------------------------------------------     Agregacja     ---------------------------------------------------------
	//---------------------------------------------------------------------------------------------------------------------------------------

	k = 0;
	for (int l = 0; l < (d.nL - 1); l++)
	{
		for (int m = 0; m < (d.nH - 1); m++)
		{



			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					this->gH[this->arrOfElements[k].nodesOfElement[i]][this->arrOfElements[k].nodesOfElement[j]] += this->arrOfElements[k].macierzHOstateczna[i][j];
					this->gC[this->arrOfElements[k].nodesOfElement[i]][this->arrOfElements[k].nodesOfElement[j]] += this->arrOfElements[k].macierzC[i][j];

				}
				this->gP[this->arrOfElements[k].nodesOfElement[i]] += this->arrOfElements[k].wektorP[i];
			}

			k++;
		}
	}// --------------------------- </> agr
	
	cout << "Glob HC : \n\n";
	for (int i = 0; i < numberOfNodes; i++) {
		for (int j = 0; j < numberOfNodes; j++) {
			gHC[i][j] = gH[i][j] + (gC[i][j] / d.stepTime);
			//cout << gHC[i][j] << "   ";
		}
		//cout << endl;
	}

	double *tc = new double[numberOfNodes];

	for (int i = 0; i < numberOfNodes; i++) {
		double sum = 0;
		for (int j = 0; j < numberOfNodes; j++) {
			sum += gC[i][j] * 100;    //TODO : here! tempetarure on start
		}
		tc[i] = sum;
	}











	//TODO: JAKAS FUNKCJA WYPISUJACA
	//---- </>tutaj test dla danych z zajec:

	for (int c = 0; c < numberOfElements; c++)
	{
		//cout << "Element " << c << "**********************************************************************" << endl;

		/*
		cout << "macierze: dN/d Ksi \t i \t dN/d Eta" << endl;
		for (int i = 0; i < 4; i++) // i determinuje nr funkcji ksztaltu
		{
			cout << this->arrOfElements[c].pochodnePoKsi[i][0] << "\t";
			cout << this->arrOfElements[c].pochodnePoKsi[i][1] << "\t";
			cout << this->arrOfElements[c].pochodnePoKsi[i][2] << "\t";
			cout << this->arrOfElements[c].pochodnePoKsi[i][3] << "\t \t";

			cout << this->arrOfElements[c].pochodnePoEta[i][0] << "\t";
			cout << this->arrOfElements[c].pochodnePoEta[i][1] << "\t";
			cout << this->arrOfElements[c].pochodnePoEta[i][2] << "\t";
			cout << this->arrOfElements[c].pochodnePoEta[i][3] << "\t";
			cout << endl;
		}
		for (int f = 0; f < 4; f++)
		{
			cout << "kom J " << this->arrOfElements[c].macierzJ[f][0] << "\t" << this->arrOfElements[c].macierzJ[f][1] << "\t" << this->arrOfElements[c].macierzJ[f][2] << "\t" << this->arrOfElements[c].macierzJ[f][3] << endl << endl;
		}

		cout << endl << "det:" << endl;
		for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
		{
			cout << this->arrOfElements[c].macierzDetJ[i] << "\t";
		}

		cout << endl << " odwr j:" << endl;
		for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
		{
			cout << this->arrOfElements[c].macierzOdwrJ[i][0] << "\t";
			cout << this->arrOfElements[c].macierzOdwrJ[i][1] << "\t";
			cout << this->arrOfElements[c].macierzOdwrJ[i][2] << "\t";
			cout << this->arrOfElements[c].macierzOdwrJ[i][3] << endl;
		}

		cout << endl << " pochodne po x" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << this->arrOfElements[c].pochodnePoX[i][0] << "\t";
			cout << this->arrOfElements[c].pochodnePoX[i][1] << "\t";
			cout << this->arrOfElements[c].pochodnePoX[i][2] << "\t";
			cout << this->arrOfElements[c].pochodnePoX[i][3] << endl;
		}

		cout << endl << " pochodne po y" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << this->arrOfElements[c].pochodnePoY[i][0] << "\t";
			cout << this->arrOfElements[c].pochodnePoY[i][1] << "\t";
			cout << this->arrOfElements[c].pochodnePoY[i][2] << "\t";
			cout << this->arrOfElements[c].pochodnePoY[i][3] << endl;
		}
		//@@

		cout << endl << " tx4" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << this->arrOfElements[c].macierzTX4[i][0] << "\t";
			cout << this->arrOfElements[c].macierzTX4[i][1] << "\t";
			cout << this->arrOfElements[c].macierzTX4[i][2] << "\t";
			cout << this->arrOfElements[c].macierzTX4[i][3] << endl;
		}
		cout << endl << " ty4" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << this->arrOfElements[c].macierzTY4[i][0] << "\t";
			cout << this->arrOfElements[c].macierzTY4[i][1] << "\t";
			cout << this->arrOfElements[c].macierzTY4[i][2] << "\t";
			cout << this->arrOfElements[c].macierzTY4[i][3] << endl;
		}

		cout << endl << " podcalk 4" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << this->arrOfElements[c].macierzPodcalkowa4[i][0] << "\t";
			cout << this->arrOfElements[c].macierzPodcalkowa4[i][1] << "\t";
			cout << this->arrOfElements[c].macierzPodcalkowa4[i][2] << "\t";
			cout << this->arrOfElements[c].macierzPodcalkowa4[i][3] << endl;
		}
	
		cout << endl << " macierz H" << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << this->arrOfElements[c].macierzH[i][0] << "\t";
			cout << this->arrOfElements[c].macierzH[i][1] << "\t";
			cout << this->arrOfElements[c].macierzH[i][2] << "\t";
			cout << this->arrOfElements[c].macierzH[i][3] << endl;
		}
		

		for (int i = 0; i < 4; i++) // i - punkt calkowania, 0123 - funkcja ksztaltu
		{
			cout << this->arrOfElements[c].macierzN[i][0]  << "\t";
			cout << this->arrOfElements[c].macierzN[i][1]  << "\t";
			cout << this->arrOfElements[c].macierzN[i][2]  << "\t";
			cout << this->arrOfElements[c].macierzN[i][3]  << endl;
		}
		
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				//cout << this->arrOfElements[c].macierzNN1[i][j] << "\t";
				cout << this->arrOfElements[c].macierzNN2[i][j] << "\t";
				//cout << this->arrOfElements[c].macierzNN3[i][j] << "\t";
				//cout << this->arrOfElements[c].macierzNN4[i][j] << endl;
			}
			cout << endl;
		}
		
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << this->arrOfElements[c].macierzC[i][j] << "\t";
			}
			cout << endl;
		}

cout << this->arrOfNodes[this->arrOfElements[6].nodesOfElement[2]].isOnEdge << "uaaaaaaaaaaa" << endl;

cout << "powierzchnie elementu  "<< c << " : " ;
cout << this->arrOfElements[6].surface[0].isEdge() << "\t";
cout << this->arrOfElements[6].surface[1].isEdge() << "\t";
cout << this->arrOfElements[6].surface[2].isEdge() << "\t";
cout << this->arrOfElements[6].surface[3].isEdge() << "\t";


cout << " ...................................................................................\n";
for (int j = 0; j < 4; j++) {
	for (int i = 0; i < 2; i++) {
		cout << this->arrOfElements[c].surface[j].macierzN[i][0] << "\t";
		cout << this->arrOfElements[c].surface[j].macierzN[i][1] << "\t";
		cout << this->arrOfElements[c].surface[j].macierzN[i][2] << "\t";
		cout << this->arrOfElements[c].surface[j].macierzN[i][3] << "\t";
		cout << endl;
	}
}

for (int j = 0; j < 4; j++) {
for (int i = 0; i < 4; i++) {
for (int h = 0; h < 4; h++) {
cout << this->arrOfElements[0].surface[j].macierzNN2[i][h] << "\t";
}
cout << endl;
}
}
*/



cout << endl;

		cout << endl;
	}


	for (int j = 0; j < 4; j++) {
		for (int i = 0; i < 4; i++) {
			for (int h = 0; h < 4; h++) {
				//cout << this->arrOfElements[0].surface[j].macierzSumNN[i][h] << "\t";
			}
			//cout << endl;
		}
	}
	
	for (int i = 0; i < numberOfNodes; i++) {
		for (int j = 0; j < numberOfNodes; j++) {
			cout << gH[i][j] << "   ";
		}
		cout << endl;
	}

	cout << endl << endl << endl;

	for (int i = 0; i < numberOfNodes; i++) {
		for (int j = 0; j < numberOfNodes; j++) {
			cout << gC[i][j] << "   ";
		}
		cout << endl;
	}

	


	cout << endl << endl << endl;

	for (int i = 0; i < numberOfNodes; i++) {
		gPC[i] = gP[i] - (tc[i] / d.stepTime);
		cout << gP[i] - (tc[i]/d.stepTime)<< "   ";

	}

	//przerzucenie na drug¹ stronê:
	for (int i = 0; i < numberOfNodes; i++) {
		gPC[i] *= (-1.0);
	}

	double *temperatures = new double[numberOfNodes];
	Gauss g;
	temperatures = g.solve(numberOfNodes, gHC, gPC);
	for (int i = 0; i < numberOfNodes; i++) {
		this ->arrOfNodes[i].t = temperatures[i];
	}


	cout << endl << endl;
	for (int i = 0; i < numberOfNodes; i++) {
		cout << "Temp dla wezla  " << i << ":   " << this->arrOfNodes[i].t << endl;
	}

	for (int i = 0; i < d.nH; i++) {
		cout << this->arrOfNodes[d.nH - 1 - i + (d.nH * 0)].t << "\t|   " << this->arrOfNodes[d.nH - 1 - i + (d.nH * 1)].t << "\t|   " << this->arrOfNodes[d.nH - 1 - i + (d.nH * 2)].t << "\t|   " << this->arrOfNodes[d.nH - 1 - i + (d.nH * 3)].t << "\t|   " << endl;
	}


	cout << endl << endl;

	double min = this->arrOfNodes[0].t;
	double max = this->arrOfNodes[0].t;
	for (int i = 1; i < numberOfNodes; i++) {
		if (min > this->arrOfNodes[i].t) {
			min = this->arrOfNodes[i].t;
		}
		if (max < this->arrOfNodes[i].t) {
			max = this->arrOfNodes[i].t;
		}
	}

	cout << "Temp min = " << min << endl << "Temp max = " << max << endl;


	for (int j = d.stepTime; j <= d.wholeTime; j+=d.stepTime) {

		cout << "po czasie " << j << " sekund:" << endl << endl;

		temperatures = g.solve(numberOfNodes, gHC, gPC);
		for (int i = 0; i < numberOfNodes; i++) {
			this->arrOfNodes[i].t = temperatures[i];
		}

		for (int i = 0; i < d.nH; i++) {
				cout << this->arrOfNodes[d.nH - 1 - i + (d.nH * 0)].t << "\t|   " << this->arrOfNodes[d.nH - 1 - i + (d.nH * 1)].t << "\t|   " << this->arrOfNodes[d.nH - 1 - i + (d.nH * 2)].t << "\t|   " << this->arrOfNodes[d.nH - 1 - i + (d.nH * 3)].t << "\t|   " << endl;
		}
	

		for (int i = 0; i < numberOfNodes; i++) {
			double sum = 0;
			for (int j = 0; j < numberOfNodes; j++) {
				sum += gC[i][j] * this->arrOfNodes[j].t; // HERE CHANGED I TO J!!!!!!!!!!!
			}
			tc[i] = sum;
		}

		for (int i = 0; i < numberOfNodes; i++) {
			gPC[i] = - (gP[i] - (tc[i] / d.stepTime));
		}

		double min = this->arrOfNodes[0].t;
		double max = this->arrOfNodes[0].t;
		for (int i = 1; i < numberOfNodes; i++) {
			if (min > this->arrOfNodes[i].t) {
				min = this->arrOfNodes[i].t;
			}
			if (max < this->arrOfNodes[i].t) {
				max = this->arrOfNodes[i].t;
			}
		}

		cout << "Temp min = " << min << endl << "Temp max = " << max << endl;

		cout << endl << " ================================" << endl;

	}












}

Grid::Grid()
{
}


Grid::~Grid()
{
}
