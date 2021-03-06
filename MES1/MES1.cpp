// MES1.cpp : nie bedzie juz potrzebny
//
// Siatka MES - odczyt danych z pliku
//
/*
#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <string>


using namespace std;



int nh; // amount of nodes in grid
int nE; // amount of elements in grid
int numberOfElements;
int numberOfNodes;
double  K = 30.0; // conductivity

//----------------------------	data from file	--------------------------------

double H; // height of the general object (on which I put the grid)
double L; // length of the general object
int nH; // number of nodes in a row vertically
int nL; // number of nodes in a row horizontaly

//----------------------------	</> data from file	--------------------------------

struct Node
{
	double x;
	double y;
	double t = 20; // Celsius degrees

};
struct xElement 
{
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


};

struct xGrid
{
	Node *arrOfNodes = new Node[numberOfNodes];
	xElement *arrOfElements = new xElement[numberOfElements];

};


void getData();
void generateGrid();
double N1(double ksi, double eta);
double N2(double ksi, double eta);
double N3(double ksi, double eta);
double N4(double ksi, double eta);


using namespace std;

int main()
{
	cout << "\n\n\nWelcome! ------------------------ Program making FEM grid ------------------------" << endl << endl;
	getData();
	//generateGrid();
	

	system("pause");
    return 0;
}

void getData()
{
	int amountOfLines = 0;
	string nameOfVariable;
	double valueDouble;
	int valueInt;

	fstream file;
	cout << "\nType file name with extention: \n"; // the file must be in THIS ccp-file's catalog
	string nameOfFile;
	cin >> nameOfFile;

	file.open(nameOfFile, ios::in);
	if (file.good() == true)
	{
		cout << "\nFile opened successfully :) \n";

		string testLine;

		while (!file.eof())
		{
			getline(file, testLine);
			++amountOfLines;

		}

		file.clear();
		file.seekg(0, ios::beg); // setting reading-pointer to the stert of the file (because after 
		//                          the loop it's at the end)


		// ------------   -----------   ------------  to modify, reading data   -----------   ----------
		file >> nameOfFile >> nameOfFile >> valueDouble;
		H = valueDouble;

		file >> nameOfFile >> nameOfFile >> valueDouble;
		L = valueDouble;

		file >> nameOfFile >> nameOfFile >> valueInt;
		nH = valueInt;

		file >> nameOfFile >> nameOfFile >> valueInt;
		nL = valueInt;

		// ------------   -----------   ------------  </>to modify, reading data   -----------   ----------

		cout << "Got the data: \n\t H = " << H << "\n\t L = " << L << "\n\t nH = " << nH;
		cout << "\n\t nL = " << nL << endl << endl;

		generateGrid();

	}
	else cout << "\nERROR - Can't open the file! :( \n";



}

void generateGrid()
{
	numberOfElements = (nH - 1) * (nL - 1);
	numberOfNodes = nH * nL;

	struct Grid
	{
		Node *arrOfNodes = new Node[numberOfNodes];
		xElement *arrOfElements = new xElement[numberOfElements];

	};
	Grid testGrid;


	//----------- ta pętla jedzie po wszystkich elementach siatki
	int k = 0;
	for (int j = 0; j < (nL - 1); j++)
	{
		for (int i = 0; i < (nH - 1); i++)
		{

			testGrid.arrOfElements[k].nodesOfElement[0] = k + j;
			testGrid.arrOfElements[k].nodesOfElement[1] = k + nH + j;
			testGrid.arrOfElements[k].nodesOfElement[2] = k + nH + 1 + j;
			testGrid.arrOfElements[k].nodesOfElement[3] = k + 1 + j;

			k++;
		}

	}
	//----------- </> ta pętla jedzie po wszystkich elementach siatki


	int a, b, c, d;
	a = testGrid.arrOfElements[0].nodesOfElement[0];
	b = testGrid.arrOfElements[(nL - 2) * (nH - 1)].nodesOfElement[1];
	c = testGrid.arrOfElements[(nL - 1) * (nH - 1) - 1].nodesOfElement[2];
	d = testGrid.arrOfElements[nH - 2].nodesOfElement[3];
	cout << "Test of nodes's index: " << a << " - " << b << " - "  << c << " - " << d << endl << endl;
	cout << " H= " << H << " L = " << L << endl;

	double deltaX = (L /(nL -1));
	double deltaY = (H / (nH - 1));
	cout << " deltax= " << deltaX << " deltaY = " << deltaY << endl;
	k = 0;
	for (int j = 0; j < nL; j++)
	{
		for (int i = 0; i < nH; i++)
		{
			testGrid.arrOfNodes[k].x = j * deltaX;
			testGrid.arrOfNodes[k].y = i * deltaY;

			k++;

		}

	}
	double x = testGrid.arrOfNodes[nH * (nL - 1)].x;
	double y = testGrid.arrOfNodes[numberOfNodes-1].y;
	cout << "Test of nodes's coordinates: \t x = " << x << "\t y = " << y << endl << endl;

	cout << endl << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements[0].nodesOfElement[i] << "\t";
	}
	cout << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements[12].nodesOfElement[i] << "\t";
	}
	cout << endl;

	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements[17].nodesOfElement[i] << "\t";
	}
	cout << endl;

	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements[5].nodesOfElement[i] << "\t";
	}
	cout << endl;


	cout << endl << endl << endl << " ################################################################";
	cout << endl << "Calkowanie:" << endl << endl;

	//----tutaj test dla danych z zajec:
	testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[0]].x = 0;
	testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[1]].x = 0.025;
	testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[2]].x = 0.025;
	testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[3]].x = 0;

	testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[0]].y = 0;
	testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[1]].y = 0;
	testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[2]].y = 0.025;
	testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[3]].y = 0.025;


	//----------- ta pętla jedzie po wszystkich elementach siatki
	k = 0;
	for (int j = 0; j < (nL - 1); j++)
	{
		for (int i = 0; i < (nH - 1); i++)
		{

			testGrid.arrOfElements[k].eta[0] = -(1 / sqrt(3));
			testGrid.arrOfElements[k].eta[1] = -(1 / sqrt(3));
			testGrid.arrOfElements[k].eta[2] =  (1 / sqrt(3));
			testGrid.arrOfElements[k].eta[3] =  (1 / sqrt(3));

			testGrid.arrOfElements[k].ksi[0] = -(1 / sqrt(3));
			testGrid.arrOfElements[k].ksi[1] =  (1 / sqrt(3));
			testGrid.arrOfElements[k].ksi[2] =  (1 / sqrt(3));
			testGrid.arrOfElements[k].ksi[3] = -(1 / sqrt(3));

			//macierz dN1/dKsi  i  dN1/dEta   dla czterech punktow calkowania
			for (int i = 0; i < 4; i++) // i determinuje punkt calkowania
			{
				testGrid.arrOfElements->pochodnePoKsi[i][0] = (-0.25 * (1 - testGrid.arrOfElements[k].eta[i]));
				testGrid.arrOfElements->pochodnePoKsi[i][1] = ( 0.25 * (1 - testGrid.arrOfElements[k].eta[i]));
				testGrid.arrOfElements->pochodnePoKsi[i][2] = ( 0.25 * (1 + testGrid.arrOfElements[k].eta[i]));
				testGrid.arrOfElements->pochodnePoKsi[i][3] = (-0.25 * (1 + testGrid.arrOfElements[k].eta[i]));

				testGrid.arrOfElements->pochodnePoEta[i][0] = (-0.25 * (1 - testGrid.arrOfElements[k].ksi[i]));
				testGrid.arrOfElements->pochodnePoEta[i][1] = (-0.25 * (1 + testGrid.arrOfElements[k].ksi[i]));
				testGrid.arrOfElements->pochodnePoEta[i][2] = (0.25 * (1 + testGrid.arrOfElements[k].ksi[i]));
				testGrid.arrOfElements->pochodnePoEta[i][3] = (0.25 * (1 - testGrid.arrOfElements[k].ksi[i]));
			}

			//macierze J dla 4 punktow calkowania
			for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
			{
				testGrid.arrOfElements->macierzJ[i][0] = ((testGrid.arrOfElements->pochodnePoKsi[i][0] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[0]].x) + (testGrid.arrOfElements->pochodnePoKsi[i][1] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[1]].x) + (testGrid.arrOfElements->pochodnePoKsi[i][2] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[2]].x) + (testGrid.arrOfElements->pochodnePoKsi[i][3] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[3]].x));
				testGrid.arrOfElements->macierzJ[i][1] = ((testGrid.arrOfElements->pochodnePoKsi[i][0] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[0]].y) + (testGrid.arrOfElements->pochodnePoKsi[i][1] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[1]].y) + (testGrid.arrOfElements->pochodnePoKsi[i][2] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[2]].y) + (testGrid.arrOfElements->pochodnePoKsi[i][3] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[3]].y));
				testGrid.arrOfElements->macierzJ[i][2] = ((testGrid.arrOfElements->pochodnePoEta[i][0] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[0]].x) + (testGrid.arrOfElements->pochodnePoEta[i][1] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[1]].x) + (testGrid.arrOfElements->pochodnePoEta[i][2] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[2]].x) + (testGrid.arrOfElements->pochodnePoEta[i][3] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[3]].x));
				testGrid.arrOfElements->macierzJ[i][3] = ((testGrid.arrOfElements->pochodnePoEta[i][0] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[0]].y) + (testGrid.arrOfElements->pochodnePoEta[i][1] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[1]].y) + (testGrid.arrOfElements->pochodnePoEta[i][2] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[2]].y) + (testGrid.arrOfElements->pochodnePoEta[i][3] * testGrid.arrOfNodes[testGrid.arrOfElements[0].nodesOfElement[3]].y));
			}

			//macierze det_J dla 4 punktow calkowania
			for (int i = 0; i < 4; i++) 
			{
				testGrid.arrOfElements->macierzDetJ[i] = ((testGrid.arrOfElements->macierzJ[i][0] * testGrid.arrOfElements->macierzJ[i][3]) - (testGrid.arrOfElements->macierzJ[i][1] * testGrid.arrOfElements->macierzJ[i][2]));
			}
			
			//macierze odwrotna do J dla 4 punktow calkowania
			for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
			{
				testGrid.arrOfElements->macierzOdwrJ[i][0] = (testGrid.arrOfElements->macierzJ[i][3] / testGrid.arrOfElements->macierzDetJ[i]);
				testGrid.arrOfElements->macierzOdwrJ[i][1] = - (testGrid.arrOfElements->macierzJ[i][1] / testGrid.arrOfElements->macierzDetJ[i]);
				testGrid.arrOfElements->macierzOdwrJ[i][2] = (testGrid.arrOfElements->macierzJ[i][2] / testGrid.arrOfElements->macierzDetJ[i]);
				testGrid.arrOfElements->macierzOdwrJ[i][3] = (testGrid.arrOfElements->macierzJ[i][0] / testGrid.arrOfElements->macierzDetJ[i]);
			}

			//macierz pochodnych po x i po y
			for (int i = 0; i < 4; i++) 
			{
				testGrid.arrOfElements->pochodnePoX[0][i] = (testGrid.arrOfElements->macierzOdwrJ[i][0] * testGrid.arrOfElements->pochodnePoKsi[i][0]) + (testGrid.arrOfElements->macierzOdwrJ[i][1] * testGrid.arrOfElements->pochodnePoEta[i][0]);
				testGrid.arrOfElements->pochodnePoX[1][i] = (testGrid.arrOfElements->macierzOdwrJ[i][0] * testGrid.arrOfElements->pochodnePoKsi[i][1]) + (testGrid.arrOfElements->macierzOdwrJ[i][1] * testGrid.arrOfElements->pochodnePoEta[i][1]);
				testGrid.arrOfElements->pochodnePoX[2][i] = (testGrid.arrOfElements->macierzOdwrJ[i][0] * testGrid.arrOfElements->pochodnePoKsi[i][2]) + (testGrid.arrOfElements->macierzOdwrJ[i][1] * testGrid.arrOfElements->pochodnePoEta[i][2]);
				testGrid.arrOfElements->pochodnePoX[3][i] = (testGrid.arrOfElements->macierzOdwrJ[i][0] * testGrid.arrOfElements->pochodnePoKsi[i][3]) + (testGrid.arrOfElements->macierzOdwrJ[i][1] * testGrid.arrOfElements->pochodnePoEta[i][3]);

				testGrid.arrOfElements->pochodnePoY[0][i] = (testGrid.arrOfElements->macierzOdwrJ[i][2] * testGrid.arrOfElements->pochodnePoKsi[i][0]) + (testGrid.arrOfElements->macierzOdwrJ[i][3] * testGrid.arrOfElements->pochodnePoEta[i][0]);
				testGrid.arrOfElements->pochodnePoY[1][i] = (testGrid.arrOfElements->macierzOdwrJ[i][2] * testGrid.arrOfElements->pochodnePoKsi[i][1]) + (testGrid.arrOfElements->macierzOdwrJ[i][3] * testGrid.arrOfElements->pochodnePoEta[i][1]);
				testGrid.arrOfElements->pochodnePoY[2][i] = (testGrid.arrOfElements->macierzOdwrJ[i][2] * testGrid.arrOfElements->pochodnePoKsi[i][2]) + (testGrid.arrOfElements->macierzOdwrJ[i][3] * testGrid.arrOfElements->pochodnePoEta[i][2]);
				testGrid.arrOfElements->pochodnePoY[3][i] = (testGrid.arrOfElements->macierzOdwrJ[i][2] * testGrid.arrOfElements->pochodnePoKsi[i][3]) + (testGrid.arrOfElements->macierzOdwrJ[i][3] * testGrid.arrOfElements->pochodnePoEta[i][3]);
			}

			//###
			for (int i = 0; i < 4; i++) 
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzTX1[i][j] = testGrid.arrOfElements->macierzDetJ[0] * testGrid.arrOfElements->pochodnePoX[i][0] * testGrid.arrOfElements->pochodnePoX[j][0];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzTX2[i][j] = testGrid.arrOfElements->macierzDetJ[1] * testGrid.arrOfElements->pochodnePoX[i][1] * testGrid.arrOfElements->pochodnePoX[j][1];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzTX3[i][j] = testGrid.arrOfElements->macierzDetJ[2] * testGrid.arrOfElements->pochodnePoX[i][2] * testGrid.arrOfElements->pochodnePoX[j][2];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzTX4[i][j] = testGrid.arrOfElements->macierzDetJ[3] * testGrid.arrOfElements->pochodnePoX[i][3] * testGrid.arrOfElements->pochodnePoX[j][3];
				}
			}


			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzTY1[i][j] = testGrid.arrOfElements->macierzDetJ[0] * testGrid.arrOfElements->pochodnePoY[i][0] * testGrid.arrOfElements->pochodnePoY[j][0];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzTY2[i][j] = testGrid.arrOfElements->macierzDetJ[1] * testGrid.arrOfElements->pochodnePoY[i][1] * testGrid.arrOfElements->pochodnePoY[j][1];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzTY3[i][j] = testGrid.arrOfElements->macierzDetJ[2] * testGrid.arrOfElements->pochodnePoY[i][2] * testGrid.arrOfElements->pochodnePoY[j][2];
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzTY4[i][j] = testGrid.arrOfElements->macierzDetJ[3] * testGrid.arrOfElements->pochodnePoY[i][3] * testGrid.arrOfElements->pochodnePoY[j][3];
				}
			}

			//##

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzPodcalkowa1[i][j] = K * (testGrid.arrOfElements->macierzTX1[i][j] + testGrid.arrOfElements->macierzTY1[i][j]);
					testGrid.arrOfElements->macierzPodcalkowa2[i][j] = K * (testGrid.arrOfElements->macierzTX2[i][j] + testGrid.arrOfElements->macierzTY2[i][j]);
					testGrid.arrOfElements->macierzPodcalkowa3[i][j] = K * (testGrid.arrOfElements->macierzTX3[i][j] + testGrid.arrOfElements->macierzTY3[i][j]);
					testGrid.arrOfElements->macierzPodcalkowa4[i][j] = K * (testGrid.arrOfElements->macierzTX4[i][j] + testGrid.arrOfElements->macierzTY4[i][j]);

				}
			}
			
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					testGrid.arrOfElements->macierzH[i][j] = (testGrid.arrOfElements->macierzPodcalkowa1[i][j] + testGrid.arrOfElements->macierzPodcalkowa2[i][j] + testGrid.arrOfElements->macierzPodcalkowa3[i][j] + testGrid.arrOfElements->macierzPodcalkowa4[i][j]);
				}
			}

			k++;// dont touch it!!! iterator elementu
		}
	}
	//----------- </> ta pętla jedzie po wszystkich elementach siatki


	

	//TODO: JAKAS FUNKCJA WYPISUJACA
	//---- </>tutaj test dla danych z zajec:

	cout << "macierze: dN/d Ksi \t i \t dN/d Eta" << endl;
	for (int i = 0; i < 4; i++) // i determinuje nr funkcji ksztaltu
	{
		cout << testGrid.arrOfElements->pochodnePoKsi[i][0] << "\t";
		cout << testGrid.arrOfElements->pochodnePoKsi[i][1] << "\t";
		cout << testGrid.arrOfElements->pochodnePoKsi[i][2] << "\t";
		cout << testGrid.arrOfElements->pochodnePoKsi[i][3] << "\t \t";

		cout << testGrid.arrOfElements->pochodnePoEta[i][0] << "\t";
		cout << testGrid.arrOfElements->pochodnePoEta[i][1] << "\t";
		cout << testGrid.arrOfElements->pochodnePoEta[i][2] << "\t";
		cout << testGrid.arrOfElements->pochodnePoEta[i][3] << "\t";
		cout << endl;
	}
	for (int f = 0; f < 4; f++)
	{
		cout << "kom J " << testGrid.arrOfElements->macierzJ[f][0] << "\t" << testGrid.arrOfElements->macierzJ[f][1] << "\t" << testGrid.arrOfElements->macierzJ[f][2] << "\t" << testGrid.arrOfElements->macierzJ[f][3] << endl << endl;
	}

	cout << endl << "det:" << endl;
	for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
	{
		cout << testGrid.arrOfElements->macierzDetJ[i] << "\t";
	}

	cout << endl << " odwr j:" << endl;
	for (int i = 0; i < 4; i++) // i determinuje punkt calkowania;   0,1 - I rzad macierzy 2x2   2,3 - drugi rzad macierzy 2x2
	{
		cout << testGrid.arrOfElements->macierzOdwrJ[i][0] << "\t";
		cout << testGrid.arrOfElements->macierzOdwrJ[i][1] << "\t";
		cout << testGrid.arrOfElements->macierzOdwrJ[i][2] << "\t";
		cout << testGrid.arrOfElements->macierzOdwrJ[i][3] << endl;
	}

	cout << endl << " pochodne po x" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements->pochodnePoX[i][0] << "\t";
		cout << testGrid.arrOfElements->pochodnePoX[i][1] << "\t";
		cout << testGrid.arrOfElements->pochodnePoX[i][2] << "\t";
		cout << testGrid.arrOfElements->pochodnePoX[i][3] << endl;
	}

	cout << endl << " pochodne po y" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements->pochodnePoY[i][0] << "\t";
		cout << testGrid.arrOfElements->pochodnePoY[i][1] << "\t";
		cout << testGrid.arrOfElements->pochodnePoY[i][2] << "\t";
		cout << testGrid.arrOfElements->pochodnePoY[i][3] << endl;
	}
	//@@

	cout << endl << " tx4" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements->macierzTX4[i][0] << "\t";
		cout << testGrid.arrOfElements->macierzTX4[i][1] << "\t";
		cout << testGrid.arrOfElements->macierzTX4[i][2] << "\t";
		cout << testGrid.arrOfElements->macierzTX4[i][3] << endl;
	}
	cout << endl << " ty4" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements->macierzTY4[i][0] << "\t";
		cout << testGrid.arrOfElements->macierzTY4[i][1] << "\t";
		cout << testGrid.arrOfElements->macierzTY4[i][2] << "\t";
		cout << testGrid.arrOfElements->macierzTY4[i][3] << endl;
	}

	cout << endl << " podcalk 4" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements->macierzPodcalkowa4[i][0] << "\t";
		cout << testGrid.arrOfElements->macierzPodcalkowa4[i][1] << "\t";
		cout << testGrid.arrOfElements->macierzPodcalkowa4[i][2] << "\t";
		cout << testGrid.arrOfElements->macierzPodcalkowa4[i][3] << endl;
	}

	cout << endl << " macierz H" << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << testGrid.arrOfElements->macierzH[i][0] << "\t";
		cout << testGrid.arrOfElements->macierzH[i][1] << "\t";
		cout << testGrid.arrOfElements->macierzH[i][2] << "\t";
		cout << testGrid.arrOfElements->macierzH[i][3] << endl;
	}

	cout << endl;

}

double N1(double ksi, double eta)
{
	return double (0.25*(1-ksi)*(1-eta));
}
double N2(double ksi, double eta)
{
	return double(0.25*(1 + ksi)*(1 - eta));
}
double N3(double ksi, double eta)
{
	return double(0.25*(1 + ksi)*(1 + eta));
}
double N4(double ksi, double eta)
{
	return double(0.25*(1 - ksi)*(1 + eta));
}


*/