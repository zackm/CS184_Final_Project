#include "Vec3.h"

#pragma once
class Particle{
public:
	Vec3 position, velocity; //in units meters and meters/sec
	float mass; //in units kg

	Particle(void){mass = .001;};
	Particle(Vec3,Vec3,float);
	Particle(float,float,float);//temporary
	~Particle(void){}
};

