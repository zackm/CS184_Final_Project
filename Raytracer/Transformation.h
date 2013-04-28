#pragma once
#include "glm/glm.hpp"

#include <vector>

using namespace std;

class Transformation {
public:
	glm::mat4 m;
	glm::mat4 minvt;
	glm::mat4 minv;

	Transformation(){};
	Transformation(vector<glm::mat4>);

	glm::vec3 object_point(glm::vec3);
	glm::vec3 object_vector(glm::vec3);
	glm::vec3 world_point(glm::vec3);
	glm::vec3 world_vector(glm::vec3);
	glm::vec3 world_normal(glm::vec3);
};