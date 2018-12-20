#pragma once
class Data
{
public:

	double H; // height of the general object (on which I put the grid)
	double L; // length of the general object
	int nH; // number of nodes in a row vertically
	int nL; // number of nodes in a row horizontaly
	double  K = 30.0; // conductivity

	Data();
	Data(double H, double L, int nH, int nL, double K);
	~Data();
};



Data::Data()
{
	this->H = 0;
	this->L = 0;
	this->nH = 0;
	this->nL = 0;
	this->K = 30.0;
}


Data::Data(double H, double L, int nH, int nL, double K) //rozszerzac te dane kiedy trzeba
{
	this->H = H;
	this->L = L;
	this->nH = nH;
	this->nL = nL;
	this->K = 30.0;

}


Data::~Data()
{
}
