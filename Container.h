
#include "Vec3.h"

#include "Particle.h"

#pragma once
class Container{
	//assume rectangular and axis aligned
public:
	Vec3 max,min; //corners
	
	bool in_container(Particle*,float); //returns true if particle in container, else it reflects particle and returns false.
	Container(Vec3,Vec3);
	Container(void){};
	~Container(void){};
};

