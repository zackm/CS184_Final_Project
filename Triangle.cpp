#include "Triangle.h"

Triangle::Triangle(Vec3 alpha, Vec3 beta, Vec3 gamma){
	a = alpha;
	b = beta;
	c = gamma;
}

Triangle::Triangle(Vec3 vec1, Vec3 vec2, Vec3 vec3, float w1, float w2, float w3, float tol){
	//linearly interpolates the vectors to make a triangle, uses tol to determine lenth.
	// we are essentially solving for t in the equation (1-t)*vec_i + t*vec_j = tol

	a = vec1;

	if(w2>=tol){
		b = vec2;
	}else{
		b = vec1 + (vec2-vec1)*(tol-w1)/(w2-w1);
	}

	if(w3>=tol){
		c = vec3;
	}else{
		c = vec1 + (vec3-vec1)*(tol-w1)/(w3-w1);
	}
}