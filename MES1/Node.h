#pragma once
class Node  //bool czy jest brzegiem
{
public:

	double x;
	double y;
	bool isOnEdge;
	double t; // Celsius degrees

	Node();
	~Node();
};



Node::Node()
{
	this->t = 100;
	this->isOnEdge = false;
}


Node::~Node()
{
}
