#include "Shape.h"
#include "RParticle.h"
#include "glm/glm.hpp"

#pragma once
class ParticleBlob : public Shape{
public:
	//Three vertices denote the triangle
	vector<RParticle*> particles;
	float t_step;
	float tolerance;
	float boundary_density;

	float kernel(glm::vec3,glm::vec3);//density kernel for ray marching
	glm::vec3 default_gradient(glm::vec3,glm::vec3);
	float density_at_point(glm::vec3);
	glm::vec3 normal_at_point(glm::vec3 position);
	glm::vec3 normal(glm::vec3);
	glm::vec3 get_normal(glm::vec3);

	bool intersect(Ray&, float*, LocalGeo*,bool*);
	bool intersect(Ray&);
	BRDF get_brdf();

	ParticleBlob(void){t_step = .01;tolerance = 500.0f; boundary_density = 1000.0f;};
	ParticleBlob(vector<RParticle*> arg_particles,glm::vec3 a,glm::vec3 d,glm::vec3 s,glm::vec3 r,glm::vec3 e,float sp){
		particles = arg_particles;
		t_step = .01;
		tolerance = 500.0f; 
		boundary_density = 1000.0f;
		brdf.ka = a;
		brdf.kd = d;
		brdf.ks = s;
		brdf.shiny = sp;
		brdf.kr = r;
		brdf.ke = e;
	};
};