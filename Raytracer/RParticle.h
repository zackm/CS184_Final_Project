#include "glm/glm.hpp"

#pragma once
class RParticle{
	//Class used only for ray tracing, not for SPH simulation.
public:
	glm::vec3 position;
	float mass,density;

	RParticle(void){};
	RParticle(float,float,float,float);
};