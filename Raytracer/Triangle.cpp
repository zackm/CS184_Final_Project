#include "Triangle.h"

#include "glm/glm.hpp"

#pragma once
#include <limits>

#pragma once
#include "Transformation.h"

#include <iostream>

using namespace std;

Triangle::Triangle(glm::vec3 arg_a,glm::vec3 arg_b,glm::vec3 arg_c,glm::vec3 ka,glm::vec3 d,glm::vec3 s,
				   glm::vec3 r,glm::vec3 e,float sp){
					   a = arg_a;
					   b = arg_b;
					   c = arg_c;
					   brdf.ka = ka;
					   brdf.kd = d;
					   brdf.ks = s;
					   brdf.kr = r;
					   brdf.ke = e;
					   brdf.shiny = sp;
					   trinormal = false;
}

Triangle::Triangle(glm::vec3 arg_a,glm::vec3 arg_b,glm::vec3 arg_c,glm::vec3 ka,glm::vec3 d,glm::vec3 s,
				   glm::vec3 r,glm::vec3 e,float sp,glm::vec3 a_norm_arg, glm::vec3 b_norm_arg, glm::vec3 c_norm_arg){
					  
					   trinormal = true;
					   a = arg_a;
					   b = arg_b;
					   c = arg_c;
					   brdf.ka = ka;
					   brdf.kd = d;
					   brdf.ks = s;
					   brdf.kr = r;
					   brdf.ke = e;
					   brdf.shiny = sp;
					   a_norm = a_norm_arg;
					   b_norm = b_norm_arg;
					   c_norm = c_norm_arg;
}

BRDF Triangle::get_brdf() {
	return brdf;
}

bool Triangle::intersect(Ray& ray, float* thit, LocalGeo* local){
	glm::vec3 vec1 = b-a;
	glm::vec3 vec2 = c-a;
	glm::vec3 vec3 = -ray.direction;

	float alpha, beta, gamma;

	glm::vec3 matrix_point = ray.position-a;

	glm::mat3 bary_matrix(vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2],vec3[0],vec3[1],vec3[2]);

	float determ = glm::determinant(bary_matrix);
	if (determ==0){
		return false;
	}
	glm::mat3 matrix_inv = glm::inverse(bary_matrix);

	glm::vec3 sol = matrix_inv*matrix_point;
	beta = sol[0];
	gamma = sol[1];
	*thit = sol[2];
	alpha = 1-beta-gamma;

	if (*thit<ray.t_min || *thit>ray.t_max){
		return false;
	}

	if (beta>=0 && gamma>=0 && (beta+gamma<=1)){ //<= is to ensure alpha is also >=0.
		glm::vec3 normal;
		if (trinormal){
			normal = alpha*a_norm+beta*b_norm+gamma*c_norm;
		}else{
			normal = glm::cross(vec1,vec2);
		}

		float mag = glm::dot(normal,normal);
		if (mag>0){
			normal /= glm::sqrt(mag);
		}

		glm::vec3 point = ray.position+(*thit)*ray.direction;

		LocalGeo temp_local(point, normal);

		*local = temp_local;
		return true;
	}else{
		return false;
	}
}

bool Triangle::intersect(Ray& ray){
	float thit = 0;
	glm::vec3 vec1 = b-a;
	glm::vec3 vec2 = c-a;
	glm::vec3 vec3 = -ray.direction;

	float beta, gamma;

	glm::vec3 matrix_point = ray.position-a;

	glm::mat3 bary_matrix(vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2],vec3[0],vec3[1],vec3[2]);

	float determ = glm::determinant(bary_matrix);
	if (determ==0){
		return false;
	}
	glm::mat3 matrix_inv = glm::inverse(bary_matrix);

	glm::vec3 sol = matrix_inv*matrix_point;

	beta = sol[0];
	gamma = sol[1];
	thit = sol[2];

	if (thit<ray.t_min || thit>ray.t_max){
		return false;
	}

	if (beta>=0 && gamma>=0 && (beta+gamma<=1)){ //<= is to ensure alpha is also >=0.
		return true;
	}else{
		return false;
	}
}