#pragma once
#include "Light.h"

class PointLight : public Light{
public:
	void generateLightRay(LocalGeo& local,Ray* lray,glm::vec3* lcolor);
	PointLight(void){t_min = 0;t_max = 100;};
	PointLight(glm::vec3,glm::vec3, Transformation);
	~PointLight(void){};
};