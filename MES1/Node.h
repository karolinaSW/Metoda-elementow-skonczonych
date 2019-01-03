#pragma once
class Node  //bool czy jest brzegiem
{
public:

	double x;
	double y;
	bool isOnEdge;
	double t = 20; // Celsius degrees

	Node();
	~Node();
};



Node::Node()
{
	this->isOnEdge = false;
}


Node::~Node()
{
}
