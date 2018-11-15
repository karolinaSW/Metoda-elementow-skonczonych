// MES1.cpp : Defines the entry point for the console application.
//
// Siatka MES - odczyt danych z pliku

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <string>


using namespace std;



int nh; // amount of nodes in grid
int nE; // amount of elements in grid
int numberOfElements;
int numberOfNodes;

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
struct Element 
{
	int nodesOfElement[4]; // 4-element array of nodes that make an element; inside are numbers of index of nodes
};


void getData();
void generateGrid();

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
		Element *arrOfElements = new Element[numberOfElements];

	};
	Grid testGrid;

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

}
