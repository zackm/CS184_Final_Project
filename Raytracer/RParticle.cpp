#include "RParticle.h"

RParticle::RParticle(float x,float y,float z,float mass_arg){
	position = glm::vec3(x,y,z);
	mass = mass_arg;
}