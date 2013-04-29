#pragma once
#include "glm/glm.hpp"

class Ray {
public:
	glm::vec3 position, direction;
	float t_min, t_max, index_of_refraction;
	bool inside_shape;

	Ray(){};
	Ray(glm::vec3, glm::vec3, float, float,float);
};