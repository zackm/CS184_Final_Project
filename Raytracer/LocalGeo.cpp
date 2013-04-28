#include "LocalGeo.h"
#include "glm/glm.hpp"

LocalGeo::LocalGeo(glm::vec3 pt,glm::vec3 nm){
	point = pt;

	float norm = glm::sqrt(glm::dot(nm,nm));

	if (norm>0){
		nm /= norm;
	}
	normal = nm;
}
