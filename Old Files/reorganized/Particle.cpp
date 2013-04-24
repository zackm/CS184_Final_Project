#include "Particle.h"

Particle::Particle(float x, float y, float z){
	 Vec3 temp(x,y,z);
	 position = temp;
}

Particle::Particle(Vec3 pos,Vec3 vel, float arg_mass){
	position = pos;
	velocity = vel;
	mass = arg_mass;
}
