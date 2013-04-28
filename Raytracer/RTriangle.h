#pragma once
#include "Shape.h"

#pragma once
#include "Transformation.h"

class RTriangle : public Shape{
public:
	//Three vertices denote the triangle
	glm::vec3 a,b,c;
	glm::vec3 a_norm, b_norm, c_norm;
	bool trinormal; //true if we have normals. False otherwise.

	RTriangle(glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,float);
	RTriangle(glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,float,
			 glm::vec3,glm::vec3,glm::vec3);
	bool intersect(Ray&, float*, LocalGeo*,bool*);
	bool intersect(Ray&);
	RTriangle(void){};
	~RTriangle(void){};
	BRDF get_brdf();
};