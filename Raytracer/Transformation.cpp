#include "Transformation.h"

#include "glm/glm.hpp"

Transformation::Transformation(vector<glm::mat4> mat_vec){
	//multiply all matrices in the vector stack to make the transformation object.
	glm::mat4 obj_to_world(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
	int n = mat_vec.size();
	for (int i = 0; i<n ; i++){
		obj_to_world = obj_to_world * mat_vec[i];
	}
	m = obj_to_world;
	minv = glm::inverse(m);
	minvt = glm::transpose(minv);
}

/**************
To object space
***************/
glm::vec3 Transformation::object_point(glm::vec3 point){
	glm::vec4 pos;
	pos[0] = point[0];
	pos[1] = point[1];
	pos[2] = point[2];
	pos[3] = 1;

	pos = minv*pos;

	glm::vec3 pos_out;
	pos_out[0] = pos[0];
	pos_out[1] = pos[1];
	pos_out[2] = pos[2];

	return pos_out;
}

glm::vec3 Transformation::object_vector(glm::vec3 vector){
	glm::vec4 dir;

	dir[0] = vector[0];
	dir[1] = vector[1];
	dir[2] = vector[2];
	dir[3] = 0;

	dir = minv*dir;

	glm::vec3 dir_out;
	dir_out[0] = dir[0];
	dir_out[1] = dir[1];
	dir_out[2] = dir[2];

	return dir_out;
}

/*************
To world space
**************/
glm::vec3 Transformation::world_point(glm::vec3 point ){
	glm::vec4 pos;
	pos[0] = point[0];
	pos[1] = point[1];
	pos[2] = point[2];
	pos[3] = 1;

	pos = m*pos;

	glm::vec3 pos_out;
	pos_out[0] = pos[0];
	pos_out[1] = pos[1];
	pos_out[2] = pos[2];

	return pos_out;
}

glm::vec3 Transformation::world_vector(glm::vec3 vector){
	glm::vec4 dir;

	dir[0] = vector[0];
	dir[1] = vector[1];
	dir[2] = vector[2];
	dir[3] = 0;

	dir = m*dir;

	glm::vec3 dir_out;
	dir_out[0] = dir[0];
	dir_out[1] = dir[1];
	dir_out[2] = dir[2];

	return dir_out;
}

glm::vec3 Transformation::world_normal(glm::vec3 vector){
	glm::vec4 dir;

	dir[0] = vector[0];
	dir[1] = vector[1];
	dir[2] = vector[2];
	dir[3] = 0;

	dir = minvt*dir;

	glm::vec3 dir_out;
	dir_out[0] = dir[0];
	dir_out[1] = dir[1];
	dir_out[2] = dir[2];

	float norm = glm::dot(dir_out,dir_out);
	if (norm>0){
		dir_out /= glm::sqrt(norm);
	}

	return dir_out;
}