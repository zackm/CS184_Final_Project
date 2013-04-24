#include "Vec3.h"

#pragma once
class Triangle{
public:
	Vec3 a,b,c,normal;

	Triangle(Vec3,Vec3,Vec3);
	Triangle(Vec3,Vec3,Vec3,float,float,float,float);
	Triangle(void){};
	~Triangle(void){};
};

