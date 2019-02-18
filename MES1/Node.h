#pragma once
class Node  
{
public:

	double x;
	double y;
	bool isOnEdge; //bool - if is on edge
	double t; // Celsius degrees

	Node();
	~Node();
};

Node::Node()
{
	this->t = 100;  //temperature of heated object on the begining  
	this->isOnEdge = false;
}

Node::~Node()
{
}
