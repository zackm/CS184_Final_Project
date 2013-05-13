#include "Particle.h"

/*
Used only for SPH simulation
*/
Particle::Particle(float x, float y, float z){
	 Vec3 temp(x,y,z);
	 position = temp;
}

Particle::Particle(Vec3 pos,Vec3 vel, float arg_mass,float density_arg){
	position = pos;
	velocity = vel;
	mass = arg_mass;
	density = density_arg;
}
