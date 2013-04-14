
#include "Vec3.h"

#include "Particle.h"

#pragma once
class Container{
	//assume rectangular and axis aligned
public:
	Vec3 max,min;
	
	bool in_container(Particle*);
	Container(Vec3,Vec3);
	Container(void);
	~Container(void);
};

