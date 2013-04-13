#include "Particle.h"

Particle::Particle(float x, float y, float z){
	 Vec3 temp(x,y,z);
	 position = temp;
}

Particle::Particle(void)
{
}


Particle::~Particle(void)
{
}
