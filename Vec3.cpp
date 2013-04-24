#include "Vec3.h"

#include <cmath>

Vec3::Vec3(float arg_x,float arg_y,float arg_z){
	x = arg_x;
	y = arg_y;
	z = arg_z;
}

Vec3 Vec3::operator + (Vec3 arg_vec){
	Vec3 new_vec;
	new_vec.x = x+arg_vec.x;
	new_vec.y = y+arg_vec.y;
	new_vec.z = z+arg_vec.z;
	return new_vec;
}

Vec3 Vec3::operator - (Vec3 arg_vec){
	Vec3 new_vec;
	new_vec.x = x-arg_vec.x;
	new_vec.y = y-arg_vec.y;
	new_vec.z = z-arg_vec.z;
	return new_vec;
}

void Vec3::operator += (Vec3 arg_vec){
	x += arg_vec.x;
	y += arg_vec.y;
	z += arg_vec.z;
}

Vec3 Vec3::operator *(float t){
	Vec3 new_vec;
	new_vec.x = x*t;
	new_vec.y = y*t;
	new_vec.z = z*t;
	return new_vec;
}

Vec3 Vec3::operator /(float t){
	Vec3 new_vec;
	new_vec.x = x/t;
	new_vec.y = y/t;
	new_vec.z = z/t;
	return new_vec;
}

/*
simple dot product between two vectors.
*/
float dot(Vec3 v1, Vec3 v2){
	return (v1.x*v2.x)+(v1.y*v2.y)+(v1.z*v2.z);
}

Vec3 cross(Vec3 first,Vec3 second){
	Vec3 new_vec;

	new_vec.x = first.y*second.z-first.z*second.y;
	new_vec.y = first.z*second.x-first.x*second.z;
	new_vec.z = first.x*second.y-first.y*second.x;

	float mag = dot(new_vec,new_vec);

	if(mag>0){
		new_vec = new_vec / sqrt(mag);
	}

	return new_vec;
}