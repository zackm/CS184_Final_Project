#include "glm/glm.hpp"

#pragma once
class RParticle{
public:
	glm::vec3 position;
	float mass,density;

	RParticle(void){};
	RParticle(float,float,float,float);
};