#include "BRDF.h"

BRDF::BRDF(glm::vec3 arg_ka,glm::vec3 arg_kd,glm::vec3 arg_ks,glm::vec3 arg_kr,glm::vec3 arg_ke,float s){
	ka = arg_ka;
	kd = arg_kd;
	ks = arg_ks;
	kr = arg_kr;
	ke = arg_ke;
	shiny = s;
}

BRDF::BRDF(void){
	glm::vec3 v(0,0,0);
	ka = kd = ks = kr = ke = v;
	shiny = 50;
}