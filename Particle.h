#include "Vec3.h"

#pragma once
#include <vector>

#pragma once
class Particle{
public:
	Vec3 position, velocity, acceleration; //in units meters and meters/sec
	float mass, density, pressure, thermal_energy; //in units kg
    int box; // index of containing box (for neighbor function)
    std::vector<int> neighbors; // vector of neighboring particle numbers (for neighbor function)

	Particle(void){mass = .001;};
	Particle(Vec3,Vec3,float);
	Particle(float,float,float);//temporary
	~Particle(void){}
    int num_neighbors(void) { return neighbors.size(); }; // returns the number of neighbors
};

