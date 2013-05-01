#include "ParticleBlob.h"
#include "../Vec3.h"
#include "../Particle.h"

const float PI = 3.1415926;
const float H = .0625;
const float SUPPORT_RADIUS = 2*H;
const float EPS = .001;
const Vec3 MAX(.5,.5,.5);
const Vec3 MIN(0,0,0);
const bool RAY_TRACING = true;

//ParticleBlob::ParticleBlob(){
//}

float ParticleBlob::kernel(glm::vec3 r_i,glm::vec3 r_j){
	glm::vec3 diff_vec = r_i-r_j;
	float mag = glm::dot(diff_vec,diff_vec);

	return (315/(64*PI*glm::pow(H,9.0f)))*glm::pow((H*H - mag),3.0f);
}

glm::vec3 ParticleBlob::default_gradient(glm::vec3 r_i,glm::vec3 r_j){
	glm::vec3 diff_vec = r_i-r_j;
	float mag = glm::dot(diff_vec,diff_vec);

	float coeff = -945.0f/(32.0f*PI*pow(H,9.0f));

	return diff_vec * coeff * pow((H*H - mag),2.0f);
}

float ParticleBlob::density_at_point(glm::vec3 position){
	float dense_value = 0.0;
	Particle* temp_part;
	//evaluate density at point i.
	int box_num = neighbors.compute_box_num(Vec3(position.x,position.y,position.z),SUPPORT_RADIUS, MAX.x,MIN.x,RAY_TRACING);//such bad code....it makes me cry
	if(box_num>=0){
		vector<int> neighbor_vec = neighbors.box_particles[box_num];
		//evaluate density at point i.
		for(int j = 0; j<neighbor_vec.size(); j++){
			temp_part = particles[neighbor_vec[j]];
			glm::vec3 part_pos(temp_part->position.x,temp_part->position.y,temp_part->position.z);
			glm::vec3 diff_vec = position - part_pos;
			float mag = glm::dot(diff_vec,diff_vec);
			if(mag<H*H){
				float weight = kernel(position,part_pos);
				dense_value += temp_part->mass * weight;
			}

		}
	}
	return dense_value;
}

glm::vec3 ParticleBlob::normal_at_point(glm::vec3 position){
	glm::vec3 normal(0,0,0);
	Particle* temp_part;
	//evaluate density at point i.
	int box_num = neighbors.compute_box_num(Vec3(position.x,position.y,position.z),SUPPORT_RADIUS, MAX.x,MIN.x,RAY_TRACING);//such bad code....it makes me cry
	if(box_num>=0){
		vector<int> neighbor_vec = neighbors.box_particles[box_num];
		//evaluate normal at point i.
		for(int j = 0; j<neighbor_vec.size(); j++){
			temp_part = particles[neighbor_vec[j]];
			glm::vec3 part_pos(temp_part->position.x,temp_part->position.y,temp_part->position.z);
			glm::vec3 diff_vec = position - part_pos;
			float mag = glm::dot(diff_vec,diff_vec);
			if(mag<H*H){
				glm::vec3 weight = default_gradient(position,part_pos);
				normal += temp_part->mass * weight / temp_part->density;
			}
		}

		float mag = glm::dot(normal,normal);
		if(mag>0){
			normal /= glm::sqrt(mag);
		}
	}
	return -1.0f*normal;
}

bool ParticleBlob::intersect(Ray& ray, float* thit, LocalGeo* local,bool* in_shape){
	float t_start = ray.t_min;
	float t_end = glm::min(ray.t_max,5.0f);

	bool hit = false;
	glm::vec3 position;
	float dense_value = 0.0;
	Particle* temp_part;
	for(float i = t_start; i<t_end; i+= t_step){
		float dense_value = 0.0;
		//get the point
		position = ray.position + i*ray.direction;
		int box_num = neighbors.compute_box_num(Vec3(position.x,position.y,position.z),SUPPORT_RADIUS, MAX.x,MIN.x,RAY_TRACING);//such bad code....it makes me cry
		if(box_num>=0){
			vector<int> neighbor_vec = neighbors.box_particles[box_num];
			//evaluate density at point i.
			for(int j = 0; j<neighbor_vec.size(); j++){
				temp_part = particles[neighbor_vec[j]];
				glm::vec3 part_pos(temp_part->position.x,temp_part->position.y,temp_part->position.z);
				glm::vec3 diff_vec = position - part_pos;//painful to type this bad code.
				float mag = glm::dot(diff_vec,diff_vec);
				if(mag<H*H){
					float weight = kernel(position,part_pos);
					dense_value += temp_part->mass * weight;
				}
			}
		}
		//boundary condition
		if((*in_shape)){
			if(dense_value<500){
			hit = true;
			(*thit) = i;
			//interpolate position to get more accurate value.

			local->point = position;

			//compute normal
			local->normal = normal_at_point(position);;
			return hit;
			}
		}else{
			if(dense_value>500){
				hit = true;
			(*thit) = i;
			//interpolate position to get more accurate value.

			local->point = position;

			//compute normal
			local->normal = normal_at_point(position);;
			return hit;
			}
		}
	}
	return hit;
}

bool ParticleBlob::intersect(Ray& ray,bool* in_shape){
	float t_start = ray.t_min;
	float t_end = glm::min(ray.t_max,5.0f);

	bool hit = false;
	glm::vec3 position;
	float dense_value = 0.0;
	Particle* temp_part;
	for(float i = t_start; i<t_end; i+= t_step){
		float dense_value = 0.0;
		//get the point
		position = ray.position + i*ray.direction;
		int box_num = neighbors.compute_box_num(Vec3(position.x,position.y,position.z),SUPPORT_RADIUS, MAX.x,MIN.x,RAY_TRACING);//such bad code....it makes me cry
		if(box_num>=0){
			vector<int> neighbor_vec = neighbors.box_particles[box_num];
			//evaluate density at point i.
			for(int j = 0; j<neighbor_vec.size(); j++){
				temp_part = particles[neighbor_vec[j]];
				glm::vec3 part_pos(temp_part->position.x,temp_part->position.y,temp_part->position.z);
				glm::vec3 diff_vec = position - part_pos;//painful to type this bad code.
				float mag = glm::dot(diff_vec,diff_vec);
				if(mag<H*H){
					float weight = kernel(position,part_pos);
					dense_value += temp_part->mass * weight;
				}
			}
		}

		if((*in_shape)){
			if(dense_value<500){
			hit = true;
			return hit;
			}
		}else{
			if(dense_value>500){
			hit = true;
			return hit;
			}
		}
	}
	return hit;
};

BRDF ParticleBlob::get_brdf(){
	return brdf;
}