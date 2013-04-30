#include "Shape.h"
#include "RParticle.h"
#include "glm/glm.hpp"
#include "../Neighbor.h"
#include "../Particle.h"

#pragma once
class ParticleBlob : public Shape{
public:
	//Three vertices denote the triangle
	vector<Particle*> particles;
	float t_step;
	float tolerance;
	float boundary_density;
	Neighbor neighbors;

	float kernel(glm::vec3,glm::vec3);//density kernel for ray marching
	glm::vec3 default_gradient(glm::vec3,glm::vec3);
	float density_at_point(glm::vec3);
	glm::vec3 normal_at_point(glm::vec3 position);
	glm::vec3 normal(glm::vec3);
	glm::vec3 get_normal(glm::vec3);

	bool intersect(Ray&, float*, LocalGeo*,bool*);
	bool intersect(Ray&);
	BRDF get_brdf();

	ParticleBlob(void){t_step = .01;tolerance = 200.0f; boundary_density = 1000.0f;};
	ParticleBlob(vector<Particle*> arg_particles,Neighbor arg_neighbors,glm::vec3 a,glm::vec3 d,glm::vec3 s,glm::vec3 r,glm::vec3 e,float sp){
		particles = arg_particles;
		t_step = .01;
		tolerance = 200.0f; 
		boundary_density = 1000.0f;
		brdf.ka = a;
		brdf.kd = d;
		brdf.ks = s;
		brdf.shiny = sp;
		brdf.kr = r;
		brdf.ke = e;
		neighbors = arg_neighbors;
	};
};