#include "Vec3.h"

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