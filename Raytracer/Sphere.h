#include "glm/glm.hpp"

#pragma once
#include "Ray.h"

#pragma once
#include "LocalGeo.h"
#include "Shape.h"

#pragma once
#include "Transformation.h"

class Sphere: public Shape{
public:
	glm::vec3 center; //in object coordinates
	float radius;
	Transformation trans; //the transformation applied from input file

	Sphere(void);
	Sphere(glm::vec3,float,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,float,Transformation);
	~Sphere(void){};
	bool intersect(Ray& ray, float* thit, LocalGeo* local,bool*);
	bool intersect(Ray& ray,bool*);
	BRDF get_brdf();
};

