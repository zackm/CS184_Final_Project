#include "Vec3.h"

#pragma once
class Particle{
public:
	Vec3 position, velocity, pressure_gradient;
	float density, number_density, temperature, mass, pressure;

	Particle(void);
	~Particle(void);
};

