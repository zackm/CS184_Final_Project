#pragma once
#include "glm/glm.hpp"
/*
This is essentially the Material class from ray tracer design note.
We store BRDF coefficients for shading.
*/
class BRDF{
public:
	glm::vec3 ka,kd,ks,kr,ke;
	float shiny;

	BRDF(void);
	BRDF(glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,float);
	~BRDF(void){};
};