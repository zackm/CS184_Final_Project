#include "Vec3.h"

#pragma once
class Particle{
public:
	Vec3 position;
	Vec3 velocity;
	float density;
	float temperature;

	float dot(Vec3){};

	Particle(void);
	~Particle(void);
};

