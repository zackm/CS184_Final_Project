#include "PointLight.h"


void PointLight::generateLightRay(LocalGeo& local,Ray* lray,glm::vec3* lcolor){
	glm::vec3 pos = local.point;
	glm::vec3 dir = position - local.point; //position is light position.

	float dist = glm::sqrt(glm::dot(dir,dir));
	if (dist>0){
		dir /= dist; //now it is unit speed.
	}

	Ray temp_ray(pos,dir,t_min,dist,1.0f); //we use distance because ray is unit speed.

	*lray = temp_ray;

	glm::vec3 temp_color(color[0],color[1],color[2]);

	*lcolor = temp_color;
}

PointLight::PointLight(glm::vec3 pos,glm::vec3 col,Transformation tr){
	trans = tr;

	position = trans.world_point(pos);
	color = col;
	t_min = .001;
	t_max = 100000;//is set to distance when generating light ray.
}