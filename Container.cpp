#include "Container.h"

Container::Container(Vec3 max_arg,Vec3 min_arg){
	max = max_arg;
	min = min_arg;
}

bool Container::in_container(Particle *part){
	bool in_cont = true;
	float friction = 1.0f;//.85f;

	Vec3 pos = part->position;
	Vec3 vel = part->velocity;
	float mass = part->mass;

	//reflect x direction
	if (pos.x>max.x){
		pos.x = max.x-(pos.x-max.x);
		vel.x = -vel.x*friction;
		in_cont = false;
	}else if(pos.x<min.x){
		pos.x = min.x+(min.x-pos.x);
		vel.x = -vel.x*friction;
		in_cont = false;
	}

	//reflect y direction
	if (pos.y>max.y){
		pos.y = max.y-(pos.y-max.y);
		vel.y = -vel.y*friction;
		in_cont = false;
	}else if(pos.y<min.y){
		pos.y = min.y+(min.y-pos.y);
		vel.y = -vel.y*friction;
		in_cont = false;
	}

	//reflect z direction
	if (pos.z>max.z){
		pos.z = max.z-(pos.z-max.z);
		vel.z = -vel.z*friction;
		in_cont = false;
	}else if(pos.z<min.z){
		pos.z = min.z+(min.z-pos.z);
		vel.z = -vel.z*friction;
		in_cont = false;
	}

	Particle reflected_part(pos,vel,mass);

	(*part) = reflected_part;
	return in_cont;
}
