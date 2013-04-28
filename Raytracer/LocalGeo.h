#pragma once
#include "glm/glm.hpp"

class LocalGeo{
public:
	glm::vec3 point;
	glm::vec3 normal;

	LocalGeo(void){};
	LocalGeo(glm::vec3,glm::vec3);
	~LocalGeo(void){};
};
