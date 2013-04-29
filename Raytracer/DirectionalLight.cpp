#include "DirectionalLight.h"

#include <limits>

void DirectionalLight::generateLightRay(LocalGeo& local,Ray* lray,glm::vec3* lcolor){
	glm::vec3 pos = local.point;
	glm::vec3 dir = direction;

	float norm = glm::dot(dir,dir);
	if (norm>0){
		dir /= glm::sqrt(norm);
	}

	Ray temp_ray(pos,dir,t_min,t_max,1.0f);

	*lray = temp_ray;

	glm::vec3 temp_color(color[0],color[1],color[2]);

	*lcolor = temp_color;
}

DirectionalLight::DirectionalLight(glm::vec3 dir,glm::vec3 col, Transformation tr){
	trans = tr;

	direction = trans.world_vector(dir);
	color = col;
	t_min = .001;
	t_max = std::numeric_limits<float>::infinity();
}