#pragma once
#include "Shape.h"

#pragma once
#include "Transformation.h"

class Triangle : public Shape{
public:
	//Three vertices denote the triangle
	glm::vec3 a,b,c;
	glm::vec3 a_norm, b_norm, c_norm;
	bool trinormal; //true if we have normals. False otherwise.

	Triangle(glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,float);
	Triangle(glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,float,
			 glm::vec3,glm::vec3,glm::vec3);
	bool intersect(Ray&, float*, LocalGeo*);
	bool intersect(Ray&);
	Triangle(void){};
	~Triangle(void){};
	BRDF get_brdf();
};