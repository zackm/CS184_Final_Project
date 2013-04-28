#include "glm/glm.hpp"
#include "Ray.h"

Ray::Ray(glm::vec3 pos, glm::vec3 dir, float min, float max){
	position = pos;
	direction = dir;
	t_min = min;
	t_max = max;
}