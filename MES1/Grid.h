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

	
	int a, b, c, cc;
	a = this->arrOfElements[0].nodesOfElement[0];
	b = this->arrOfElements[(d.nL - 2) * (d.nH - 1)].nodesOfElement[1];
	c = this->arrOfElements[(d.nL - 1) * (d.nH - 1) - 1].nodesOfElement[2];
	cc = this->arrOfElements[d.nH - 2].nodesOfElement[3];
	cout << "Test of nodes's index: " << a << " - " << b << " - " << c << " - " << cc << endl << endl;
	cout << " H= " << d.H << " L = " << d.L << endl;
	
	int e;
	a = this->arrOfElements[0].nodesOfElement[0];
	b = this->arrOfElements[(d.nL - 2) * (d.nH - 1)].nodesOfElement[1];
	c = this->arrOfElements[(d.nL - 1) * (d.nH - 1) - 1].nodesOfElement[2];
	e = this->arrOfElements[d.nH - 2].nodesOfElement[3];
	cout << "Test of nodes's index: " << a << " - " << b << " - " << c << " - " << e << endl << endl;
	cout << " H= " << d.H << " L = " << d.L << endl;

	double deltaX = (d.L / (d.nL - 1));
	double deltaY = (d.H / (d.nH - 1));
	cout << " deltax= " << deltaX << " deltaY = " << deltaY << endl;
	k = 0;
	for (int j = 0; j < d.nL; j++)
	{
		for (int i = 0; i < d.nH; i++)
		{
			this->arrOfNodes[k].x = j * deltaX;
			this->arrOfNodes[k].y = i * deltaY;

			k++;

		}

	}
	double x = this->arrOfNodes[d.nH * (d.nL - 1)].x;
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


	cout << endl << endl << endl << " ################################################################";
	cout << endl << "Calkowanie:" << endl << endl;

	//----tutaj test dla danych z zajec:
	this->arrOfNodes[this->arrOfElements[0].nodesOfElement[0]].x = 0.0;
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
			this->arrOfElements[k].eta[2] = (1 / sqrt(3));
			this->arrOfElements[k].eta[3] = (1 / sqrt(3));

			this->arrOfElements[k].ksi[0] = -(1 / sqrt(3));
			this->arrOfElements[k].ksi[1] = (1 / sqrt(3));
			this->arrOfElements[k].ksi[2] = (1 / sqrt(3));
			this->arrOfElements[k].ksi[3] = -(1 / sqrt(3));

			//macierz dN1/dKsi  i  dN1/dEta   dla czterech punktow calkowania
			for (int i = 0; i < 4; i++) // i determinuje punkt calkowania
			{
				this->arrOfElements->pochodnePoKsi[i][0] = (-0.25 * (1 - this->arrOfElements[k].eta[i]));
				this->arrOfElements->pochodnePoKsi[i][1] = (0.25 * (1 - this->arrOfElements[k].eta[i]));
				this->arrOfElements->pochodnePoKsi[i][2] = (0.25 * (1 + this->arrOfElements[k].eta[i]));
				this->arrOfElements->pochodnePoKsi[i][3] = (-0.25 * (1 + this->arrOfElements[k].eta[i]));

				this->arrOfElements->pochodnePoEta[i][0] = (-0.25 * (1 - this->arrOfElements[k].ksi[i]));
				this->arrOfElements->pochodnePoEta[i][1] = (-0.25 * (1 + this->arrOfElements[k].ksi[i]));
				this->arrOfElements->pochodnePoEta[i][2] = (0.25 * (1 + this->arrOfElements[k].ksi[i]));
				this->arrOfElements->pochodnePoEta[i][3] = (0.25 * (1 - this->arrOfElements[k].ksi[i]));
			}

			//macierze J dla 4 punktow calkowania
			for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
			{
				this->arrOfElements->macierzJ[i][0] = ((this->arrOfElements->pochodnePoKsi[i][0] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[0]].x) + (this->arrOfElements->pochodnePoKsi[i][1] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[1]].x) + (this->arrOfElements->pochodnePoKsi[i][2] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[2]].x) + (this->arrOfElements->pochodnePoKsi[i][3] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[3]].x));
				this->arrOfElements->macierzJ[i][1] = ((this->arrOfElements->pochodnePoKsi[i][0] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[0]].y) + (this->arrOfElements->pochodnePoKsi[i][1] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[1]].y) + (this->arrOfElements->pochodnePoKsi[i][2] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[2]].y) + (this->arrOfElements->pochodnePoKsi[i][3] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[3]].y));
				this->arrOfElements->macierzJ[i][2] = ((this->arrOfElements->pochodnePoEta[i][0] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[0]].x) + (this->arrOfElements->pochodnePoEta[i][1] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[1]].x) + (this->arrOfElements->pochodnePoEta[i][2] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[2]].x) + (this->arrOfElements->pochodnePoEta[i][3] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[3]].x));
				this->arrOfElements->macierzJ[i][3] = ((this->arrOfElements->pochodnePoEta[i][0] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[0]].y) + (this->arrOfElements->pochodnePoEta[i][1] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[1]].y) + (this->arrOfElements->pochodnePoEta[i][2] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[2]].y) + (this->arrOfElements->pochodnePoEta[i][3] * this->arrOfNodes[this->arrOfElements[0].nodesOfElement[3]].y));
			}

			//macierze det_J dla 4 punktow calkowania
			for (int i = 0; i < 4; i++)
			{
				this->arrOfElements->macierzDetJ[i] = (this->arrOfElements->macierzJ[i][0] * this->arrOfElements->macierzJ[i][3]) - (this->arrOfElements->macierzJ[i][1] * this->arrOfElements->macierzJ[i][2]);
			}

			//macierze odwrotna do J dla 4 punktow calkowania
			for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
			{
				this->arrOfElements->macierzOdwrJ[i][0] =  (this->arrOfElements->macierzJ[i][3] / this->arrOfElements->macierzDetJ[i]);
				this->arrOfElements->macierzOdwrJ[i][1] = -(this->arrOfElements->macierzJ[i][1] / this->arrOfElements->macierzDetJ[i]);
				this->arrOfElements->macierzOdwrJ[i][2] = (this->arrOfElements->macierzJ[i][2] / this->arrOfElements->macierzDetJ[i]);
				this->arrOfElements->macierzOdwrJ[i][3] = (this->arrOfElements->macierzJ[i][0] / this->arrOfElements->macierzDetJ[i]);
			}

			//macierz pochodnych po x i po y
			for (int i = 0; i < 4; i++)
			{
				this->arrOfElements->pochodnePoX[0][i] = (this->arrOfElements->macierzOdwrJ[i][0] * this->arrOfElements->pochodnePoKsi[i][0]) + (this->arrOfElements->macierzOdwrJ[i][1] * this->arrOfElements->pochodnePoEta[i][0]);
				this->arrOfElements->pochodnePoX[1][i] = (this->arrOfElements->macierzOdwrJ[i][0] * this->arrOfElements->pochodnePoKsi[i][1]) + (this->arrOfElements->macierzOdwrJ[i][1] * this->arrOfElements->pochodnePoEta[i][1]);
				this->arrOfElements->pochodnePoX[2][i] = (this->arrOfElements->macierzOdwrJ[i][0] * this->arrOfElements->pochodnePoKsi[i][2]) + (this->arrOfElements->macierzOdwrJ[i][1] * this->arrOfElements->pochodnePoEta[i][2]);
				this->arrOfElements->pochodnePoX[3][i] = (this->arrOfElements->macierzOdwrJ[i][0] * this->arrOfElements->pochodnePoKsi[i][3]) + (this->arrOfElements->macierzOdwrJ[i][1] * this->arrOfElements->pochodnePoEta[i][3]);

				this->arrOfElements->pochodnePoY[0][i] = (this->arrOfElements->macierzOdwrJ[i][2] * this->arrOfElements->pochodnePoKsi[i][0]) + (this->arrOfElements->macierzOdwrJ[i][3] * this->arrOfElements->pochodnePoEta[i][0]);
				this->arrOfElements->pochodnePoY[1][i] = (this->arrOfElements->macierzOdwrJ[i][2] * this->arrOfElements->pochodnePoKsi[i][1]) + (this->arrOfElements->macierzOdwrJ[i][3] * this->arrOfElements->pochodnePoEta[i][1]);
				this->arrOfElements->pochodnePoY[2][i] = (this->arrOfElements->macierzOdwrJ[i][2] * this->arrOfElements->pochodnePoKsi[i][2]) + (this->arrOfElements->macierzOdwrJ[i][3] * this->arrOfElements->pochodnePoEta[i][2]);
				this->arrOfElements->pochodnePoY[3][i] = (this->arrOfElements->macierzOdwrJ[i][2] * this->arrOfElements->pochodnePoKsi[i][3]) + (this->arrOfElements->macierzOdwrJ[i][3] * this->arrOfElements->pochodnePoEta[i][3]);
			}

			//###
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzTX1[i][j] = this->arrOfElements->macierzDetJ[0] * this->arrOfElements->pochodnePoX[i][0] * this->arrOfElements->pochodnePoX[j][0];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzTX2[i][j] = this->arrOfElements->macierzDetJ[1] * this->arrOfElements->pochodnePoX[i][1] * this->arrOfElements->pochodnePoX[j][1];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzTX3[i][j] = this->arrOfElements->macierzDetJ[2] * this->arrOfElements->pochodnePoX[i][2] * this->arrOfElements->pochodnePoX[j][2];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzTX4[i][j] = this->arrOfElements->macierzDetJ[3] * this->arrOfElements->pochodnePoX[i][3] * this->arrOfElements->pochodnePoX[j][3];
				}
			}


			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzTY1[i][j] = this->arrOfElements->macierzDetJ[0] * this->arrOfElements->pochodnePoY[i][0] * this->arrOfElements->pochodnePoY[j][0];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzTY2[i][j] = this->arrOfElements->macierzDetJ[1] * this->arrOfElements->pochodnePoY[i][1] * this->arrOfElements->pochodnePoY[j][1];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzTY3[i][j] = this->arrOfElements->macierzDetJ[2] * this->arrOfElements->pochodnePoY[i][2] * this->arrOfElements->pochodnePoY[j][2];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzTY4[i][j] = this->arrOfElements->macierzDetJ[3] * this->arrOfElements->pochodnePoY[i][3] * this->arrOfElements->pochodnePoY[j][3];
				}
			}

			//##

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzPodcalkowa1[i][j] = d.K * (this->arrOfElements->macierzTX1[i][j] + this->arrOfElements->macierzTY1[i][j]);
					this->arrOfElements->macierzPodcalkowa2[i][j] = d.K * (this->arrOfElements->macierzTX2[i][j] + this->arrOfElements->macierzTY2[i][j]);
					this->arrOfElements->macierzPodcalkowa3[i][j] = d.K * (this->arrOfElements->macierzTX3[i][j] + this->arrOfElements->macierzTY3[i][j]);
					this->arrOfElements->macierzPodcalkowa4[i][j] = d.K * (this->arrOfElements->macierzTX4[i][j] + this->arrOfElements->macierzTY4[i][j]);

				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					this->arrOfElements->macierzH[i][j] = (this->arrOfElements->macierzPodcalkowa1[i][j] + this->arrOfElements->macierzPodcalkowa2[i][j] + this->arrOfElements->macierzPodcalkowa3[i][j] + this->arrOfElements->macierzPodcalkowa4[i][j]);
				}
			}

			k++;// dont touch it!!! iterator elementu
		}
	}
	//----------- </> ta pêtla jedzie po wszystkich elementach siatki



	
	//TODO: JAKAS FUNKCJA WYPISUJACA
	//---- </>tutaj test dla danych z zajec:

	cout << "macierze: dN/d Ksi \t i \t dN/d Eta" << endl;
	for (int i = 0; i < 4; i++) // i determinuje nr funkcji ksztaltu
	{
		cout << this->arrOfElements->pochodnePoKsi[i][0] << "\t";
		cout << this->arrOfElements->pochodnePoKsi[i][1] << "\t";
		cout << this->arrOfElements->pochodnePoKsi[i][2] << "\t";
		cout << this->arrOfElements->pochodnePoKsi[i][3] << "\t \t";

		cout << this->arrOfElements->pochodnePoEta[i][0] << "\t";
		cout << this->arrOfElements->pochodnePoEta[i][1] << "\t";
		cout << this->arrOfElements->pochodnePoEta[i][2] << "\t";
		cout << this->arrOfElements->pochodnePoEta[i][3] << "\t";
		cout << endl;
	}
	for (int f = 0; f < 4; f++)
	{
		cout << "kom J " << this->arrOfElements->macierzJ[f][0] << "\t" << this->arrOfElements->macierzJ[f][1] << "\t" << this->arrOfElements->macierzJ[f][2] << "\t" << this->arrOfElements->macierzJ[f][3] << endl << endl;
	}

	cout << endl << "det:" << endl;
	for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
	{
		cout << this->arrOfElements->macierzDetJ[i] << "\t";
	}

	cout << endl << " odwr j:" << endl;
	for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
	{
		cout << this->arrOfElements->macierzOdwrJ[i][0] << "\t";
		cout << this->arrOfElements->macierzOdwrJ[i][1] << "\t";
		cout << this->arrOfElements->macierzOdwrJ[i][2] << "\t";
		cout << this->arrOfElements->macierzOdwrJ[i][3] << endl;
	}

	cout << endl << " pochodne po x" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements->pochodnePoX[i][0] << "\t";
		cout << this->arrOfElements->pochodnePoX[i][1] << "\t";
		cout << this->arrOfElements->pochodnePoX[i][2] << "\t";
		cout << this->arrOfElements->pochodnePoX[i][3] << endl;
	}

	cout << endl << " pochodne po y" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements->pochodnePoY[i][0] << "\t";
		cout << this->arrOfElements->pochodnePoY[i][1] << "\t";
		cout << this->arrOfElements->pochodnePoY[i][2] << "\t";
		cout << this->arrOfElements->pochodnePoY[i][3] << endl;
	}
	//@@

	cout << endl << " tx4" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements->macierzTX4[i][0] << "\t";
		cout << this->arrOfElements->macierzTX4[i][1] << "\t";
		cout << this->arrOfElements->macierzTX4[i][2] << "\t";
		cout << this->arrOfElements->macierzTX4[i][3] << endl;
	}
	cout << endl << " ty4" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements->macierzTY4[i][0] << "\t";
		cout << this->arrOfElements->macierzTY4[i][1] << "\t";
		cout << this->arrOfElements->macierzTY4[i][2] << "\t";
		cout << this->arrOfElements->macierzTY4[i][3] << endl;
	}

	cout << endl << " podcalk 4" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements->macierzPodcalkowa4[i][0] << "\t";
		cout << this->arrOfElements->macierzPodcalkowa4[i][1] << "\t";
		cout << this->arrOfElements->macierzPodcalkowa4[i][2] << "\t";
		cout << this->arrOfElements->macierzPodcalkowa4[i][3] << endl;
	}
	
	cout << endl << " macierz H" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << this->arrOfElements->macierzH[i][0] << "\t";
		cout << this->arrOfElements->macierzH[i][1] << "\t";
		cout << this->arrOfElements->macierzH[i][2] << "\t";
		cout << this->arrOfElements->macierzH[i][3] << endl;
	}

	cout << endl;

}

Grid::Grid()
{
}


Grid::~Grid()
{
}
