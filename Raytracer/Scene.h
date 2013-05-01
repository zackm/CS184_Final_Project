#include "glm/glm.hpp"
#include "Light.h"
#include "Shape.h"

#ifndef __CAMERA_H__
#define __CAMERA_H__
#include "Camera.h"
#endif

#ifndef __FILM_H__
#define __FILM_H__
#include "Film.h"
#endif

#pragma once
#include <list>

using namespace std;

#pragma once
#include "BRDF.h"

class Scene {
public:
	glm::vec3 eye_position;
	glm::vec3 UL, UR, LL, LR;
	glm::vec3 ka, kd, ks, kr;
	int width, height, maxdepth;
	list<Light*> lights;
	list<Shape*> shapes;

	Scene(){};
	Scene(glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,int,int,int);
	void set_params(glm::vec3,glm::vec3,glm::vec3,glm::vec3,glm::vec3,int,int,int);
	void render(Camera, Film);
	void add_shape(Shape*);
	void add_light(Light*);

	Ray generateReflectionRay(LocalGeo&,Ray*);
	Ray generateRefractionRay(LocalGeo&,Ray*,float,float);
	float reflectance(LocalGeo&,Ray&,float,float);
	glm::vec3 trace(Ray &,int,bool);
	bool intersect_checker(Ray&,bool*);

	glm::vec3 shading(LocalGeo, BRDF, Ray, glm::vec3,glm::vec3);
	Ray createReflectRay(LocalGeo, Ray);
};