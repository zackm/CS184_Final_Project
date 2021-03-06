#pragma once
class Vec3{
public:
	float x,y,z;

	//overloads for usual additions and multiplication
	Vec3 operator + (Vec3);
	Vec3 operator - (Vec3);
	void operator += (Vec3);
	Vec3 operator * (float);
	Vec3 operator / (float);

	Vec3(void){};
	Vec3(float,float,float);

	~Vec3(void){};
};