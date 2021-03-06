#include "Vec3.h"

#pragma once
class Triangle{
public:
	Vec3 a,b,c,a_normal,b_normal,c_normal;

	Triangle(Vec3,Vec3,Vec3);
	Triangle(Vec3,Vec3,Vec3,Vec3,Vec3,Vec3); //constructor with normals
	Triangle(Vec3,Vec3,Vec3,float,float,float,float);
	Triangle(void){};
	~Triangle(void){};
};

extern float dot(Vec3,Vec3);

extern Vec3 cross(Vec3,Vec3);
