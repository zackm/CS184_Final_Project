#include "Container.h"

#include "Particle.h"

#pragma once
class Rectangle : public Container{
public:
	Vec3 max,min;

	bool in_container(Particle*,float);
	Rectangle(void){};
	Rectangle(Vec3,Vec3);
	~Rectangle(void){};
};

