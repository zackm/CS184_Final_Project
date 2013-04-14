
#include "Vec3.h"

#pragma once
class Container{
	//assume rectangular and axis aligned
public:
	Vec3 max,min;
	
	Container(Vec3,Vec3);
	Container(void);
	~Container(void);
};

