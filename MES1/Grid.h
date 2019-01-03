#pragma once

#include "Element.h"
#include "Data.h"

using namespace std;

class Grid
{
public:

	//int nh; // amount of nodes in grid
	//int nE; // amount of elements in grid
	int numberOfElements;
	int numberOfNodes;

	//Node *arrOfNodes = new Node[numberOfNodes];
	//Element *arrOfElements = new Element[numberOfElements];

	Node *arrOfNodes;
	Element *arrOfElements;

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

	cout << endl << endl << endl << " ################################################################";
	cout << endl << "Calkowanie:" << endl << endl;
	

	//----tutaj test dla danych z zajec:
	this->arrOfNodes[this->arrOfElements[0].nodesOfElement[0]].x = 0.0; //po konkretnym elemencie wyszukujê numery wêz³ów w tym elemencie i dajê je jako index do tablicy ogólmej wêz³ów
	this->arrOfNodes[this->arrOfElements[0].nodesOfElement[1]].x = 0.025;
	this->arrOfNodes[this->arrOfElements[0].nodesOfElement[2]].x = 0.025;
	this->arrOfNodes[this->arrOfElements[0].nodesOfElement[3]].x = 0.0;

	this->arrOfNodes[this->arrOfElements[0].nodesOfElement[0]].y = 0.0;
	this->arrOfNodes[this->arrOfElements[0].nodesOfElement[1]].y = 0.0;
	this->arrOfNodes[this->arrOfElements[0].nodesOfElement[2]].y = 0.025;
	this->arrOfNodes[this->arrOfElements[0].nodesOfElement[3]].y = 0.025;
	
	
	//----------- ta pêtla jedzie po wszystkich elementach siatki
	k = 0;
	for (int j = 0; j < (d.nL - 1); j++)
	{
		for (int i = 0; i < (d.nH - 1); i++)
		{

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









			k++;// dont touch it!!! iterator elementu
		}
	}
	//----------- </> ta pêtla jedzie po wszystkich elementach siatki



	















	//TODO: JAKAS FUNKCJA WYPISUJACA
	//---- </>tutaj test dla danych z zajec:

	for (int c = 0; c < numberOfElements; c++)
	{
		cout << "Element " << c << "**********************************************************************" << endl;

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
*/

		cout << endl;
	}

}

Grid::Grid()
{
}


Grid::~Grid()
{
}
