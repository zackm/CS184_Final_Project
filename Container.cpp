#include "Container.h"

const float COLLISOIN_DISTANCE = .01;

Container::Container(Vec3 max_arg,Vec3 min_arg){
	max = max_arg;
	min = min_arg;
}

bool Container::in_container(Particle *part,float t){
	//We should check time as well.

	bool in_cont = true;
	float friction = .5f;//.5f;

	Vec3 pos = part->position;
	Vec3 vel = part->velocity;
	float mass = part->mass;

	//reflect x position
	if (pos.x>max.x){
		pos.x = max.x-(pos.x-max.x);
		vel.x = -vel.x*friction;
		in_cont = false;
	}else if(pos.x<min.x){
		pos.x = min.x+(min.x-pos.x);
		vel.x = -vel.x*friction;
		in_cont = false;
	}

	//reflect y position
	if (pos.y>max.y){
		pos.y = max.y-(pos.y-max.y);
		vel.y = -vel.y*friction;
		in_cont = false;
	}else if(pos.y<min.y){
		pos.y = min.y+(min.y-pos.y);
		vel.y = -vel.y*friction;
		in_cont = false;
	}

	//reflect z position
	if (pos.z>max.z){
		pos.z = max.z-(pos.z-max.z);
		vel.z = -vel.z*friction;
		in_cont = false;
	}else if(pos.z<min.z){
		pos.z = min.z+(min.z-pos.z);
		vel.z = -vel.z*friction;
		in_cont = false;
	}

	part->position = pos;
	part->velocity = vel;
	return in_cont;


	
}