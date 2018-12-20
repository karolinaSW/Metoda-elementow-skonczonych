/// Do wprowadzenia nowej danej: 1. dodaj nowa zmienna w klasie Data, 2. dodaj wypisywanie w readData, 3. dodaj zmienna w readData, 
///4. dodaj odczyt w readData, 5. dodaj dane w pliku dane.txt

#include "stdafx.h"
#include "targetver.h"
#include <iostream>
#include <string>
#include <fstream>
#include "Grid.h"
#include "Data.h"

using namespace std;

void readData(Data &d);

//#########################################################################################################################


int main()
{
	std::cout << "\n\n\nWelcome!\n\n ------------------------ Program FEM ------------------------" << endl << endl;

	Data data(0.3, 0.1, 7, 4, 30.0); //H, L, nH, nL
	readData(data);
	Grid g;
	g.generateGrid(data);
	//generateGrid();


	system("pause");
	return 0;
}

//#########################################################################################################################

void readData(Data &d)
{
	//Czy da rade oczyscic z tych zmiennych - TAK

	/*
	double H; // height of the general object (on which I put the grid)
	double L; // length of the general object
	int nH; // number of nodes in a row vertically
	int nL; // number of nodes in a row horizontaly
	double K; // conductivity
	*/

	char option;
	cout << "From where to read data? Choose... \nFile txt  -  1 \t\t Code of the program (from class)  -  2\n ";
	cin >> option;

	switch (option) {
	case '1':
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
			file.seekg(0, ios::beg); // setting reading-pointer to the start of the file (because after 
									 //                          the loop it's at the end)


									 // ------------   -----------   ------------  (to modify) reading data   -----------   ----------
			file >> nameOfFile >> nameOfFile >> valueDouble;
			d.H = valueDouble;
			//d.H = H;

			file >> nameOfFile >> nameOfFile >> valueDouble;
			d.L = valueDouble;
			//d.L = L;

			file >> nameOfFile >> nameOfFile >> valueInt;
			d.nH = valueInt;
			//d.nH = nH;

			file >> nameOfFile >> nameOfFile >> valueInt;
			d.nL = valueInt;
			//d.nL = nL;

			file >> nameOfFile >> nameOfFile >> valueDouble;
			d.K = valueDouble;
			//d.K = K;
			// ------------   -----------   ------------  </>(to modify) reading data   -----------   ----------

			cout << "Got the data: \n\t H = " << d.H << "\n\t L = " << d.L << "\n\t nH = " << d.nH;
			cout << "\n\t nL = " << d.nL << "\n\t K = " << d.K << endl << endl;
		}
		else cout << "\nERROR - Can't open the file! :( \n";

		break;
	}
	case '2':
	{
		// write parameters (data) here
		cout << "Got the data: \n\t H = " << d.H << "\n\t L = " << d.L << "\n\t nH = " << d.nH;
		cout << "\n\t nL = " << d.nL << "\n\t K = " << d.K << endl << endl;

		break;
	}
	}
}