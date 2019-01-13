#pragma once
class Data
{
public:

	double H; // height of the general object (on which I put the grid), in real measurement units like [mm], [cm], [m]
	double L; // length of the general object, in real measurement units like [mm], [cm], [m]
	int nH; // number of nodes in a row vertically, amount - so counting from 1
	int nL; // number of nodes in a row horizontaly, amount - so counting from 1
	double  K; // conductivity
	double c; // pojemnosc cieplna
	double ro; // density
	double alfa; //convection alfa

	Data();
	Data(double H, double L, int nH, int nL, double K, double c, double ro, double alfa);
	~Data();
};



Data::Data()
{
	this->H = 0;
	this->L = 0;
	this->nH = 0;
	this->nL = 0;
	this->K = 0;
	this->c = 0;
	this->ro = 0;
	this->alfa = 0;
}


Data::Data(double H, double L, int nH, int nL, double K, double c, double ro, double alfa) //rozszerzac te dane kiedy trzeba
{
	this->H = H;
	this->L = L;
	this->nH = nH;
	this->nL = nL;
	this->K = K;
	this->c = c;
	this->ro = ro;
	this->alfa = alfa;

}


Data::~Data()
{
}
