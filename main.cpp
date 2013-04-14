#include "Particle.h"

#include <vector>

#include <GL/glut.h>
#include <GL/glu.h>

#include <stdlib.h>
#include <iostream>

#pragma once
#include "Vec3.h"

using namespace std;

/*******************
* GLOBAL VARIABLES *
*******************/

const float TIMESTEP = .01;
vector<Particle*> PARTICLES;
const int NUM_PARTICLES = 100;
const Vec3 GRAVITY(0,-9.8f,0);
const float IDEAL_DENSITY = 1.0; //for water
const float STIFFNESS = 100.0f; //no idea what it should be set to.

/************
* Overloads *
************/
//might need to overload division as well.
Vec3 operator * (float t, const Vec3& arg_vec){
	Vec3 new_vec;
	new_vec.x = arg_vec.x*t;
	new_vec.y = arg_vec.y*t;
	new_vec.z = arg_vec.z*t;
	return new_vec;
}

Vec3 operator * (const Vec3& arg_vec,float t){
	Vec3 new_vec;
	new_vec.x = arg_vec.x*t;
	new_vec.y = arg_vec.y*t;
	new_vec.z = arg_vec.z*t;
	return new_vec;
}

Vec3 operator / (float t, const Vec3& arg_vec){
	Vec3 new_vec;
	new_vec.x = arg_vec.x/t;
	new_vec.y = arg_vec.y/t;
	new_vec.z = arg_vec.z/t;
	return new_vec;
}

Vec3 operator / (const Vec3& arg_vec,float t){
	Vec3 new_vec;
	new_vec.x = arg_vec.x/t;
	new_vec.y = arg_vec.y/t;
	new_vec.z = arg_vec.z/t;
	return new_vec;
}

float dot(Vec3 v1, Vec3 v2){
	return v1.x*v2.x+v1.y*v2.y+v1.z+v2.z;
}

/******************
* Kernel Function *
******************/
float gaussian(Vec3 r_i, Vec3 r_j){
	//should add another parameter (max distance value)
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	return exp(-4.0f*mag);
}

Vec3 gaussian_grad(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float coeff = -8.0f*sqrt(mag)*exp(-4.0f*mag);

	Vec3 grad(diff_vec.x/sqrt(mag),diff_vec.y/sqrt(mag),diff_vec.z/sqrt(mag));

	return coeff*grad;
}

/********************
* Physics Functions *
********************/
Vec3 kinematic_polynomial(Vec3 acc, Vec3 vel, Vec3 pos,float t){
	return .5f*acc*t*t+vel*t+pos;
}

/*
Create a new list of particles from the old list. Then, throw the old list.
To do this, calculate all quanities in Navier-Stokes, then use timestep to
update particle location from old location and velocity.
*/
void update_particles(){
	vector<Particle*> new_particles;
	vector<float> density_list;
	vector<float> pressure_list;
	vector<Vec3> pressure_grad_list;

	//update using slow algorithm for now
	Particle* base_particle, *temp_particle, *new_particle;
	float number_density = 0, density = 0, pressure = 0;

	//Sets density at each point
	for (int i = 0; i<NUM_PARTICLES; i++){
		density = 0;
		base_particle = PARTICLES[i];

		for (int j = 0; j<NUM_PARTICLES; j++){
			//update density
			temp_particle = PARTICLES[j];
			density += temp_particle->mass*gaussian(base_particle->position,temp_particle->position);
		}
		density_list.push_back(density);
		pressure_list.push_back(STIFFNESS*(density-IDEAL_DENSITY));
	}

	//Sets pressure gradient at each point using densities from last loop
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 pressure_gradient(0,0,0);
		for (int j = 0; j<NUM_PARTICLES && j!=i; j++){
			temp_particle = PARTICLES[j];

			Vec3 weight = gaussian_grad(base_particle->position,temp_particle->position);
			pressure_gradient += temp_particle->mass * (pressure_list[j]/density_list[j])*weight; 
								 
		}
		pressure_grad_list.push_back(pressure_gradient);
	}

	//Create new particles from old particles and from pressure gradient.
	for (int i = 0; i<NUM_PARTICLES; i++){
		temp_particle = PARTICLES[i];
		Vec3 acceleration = -1.0f*pressure_grad_list[i]/density_list[i] + GRAVITY;
		Vec3 velocity = temp_particle->velocity;
		Vec3 position = temp_particle->position;
		Vec3 new_position = kinematic_polynomial(acceleration,velocity,position,TIMESTEP);
		Vec3 new_velocity = temp_particle->velocity + acceleration*TIMESTEP;
		float mass = temp_particle->mass;
		temp_particle = new Particle(new_position,new_velocity,mass);
		new_particles.push_back(temp_particle);
	}
	PARTICLES = new_particles;
}

void initScene(){
	//GLfloat light_position[] = { -1.0, -1.0, -1.0, 0.0 };
	//GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	//GLfloat mat_diffuse[] = { 0.5, 0.0, 0.7, 1.0 };
	//GLfloat mat_ambient[] = { 0.1, 0.1, 0.1, 1.0 };
	//GLfloat mat_shininess[] = { 20.0 };
	//glShadeModel(GL_SMOOTH);

	//glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	//glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	//glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	//glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	//glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	//set window and camera


	//create a list of random particles
	float x,y,z;
	for (int i = 0; i<NUM_PARTICLES; i++){
		x = float(rand())/(float(RAND_MAX));
		y = float(rand())/(float(RAND_MAX));
		z = 0.0f; //rand() % 400;
		Vec3 pos(x,y,z);
		Vec3 vel(.1,0,0);
		float mass = 1.0f;
		PARTICLES.push_back(new Particle(pos,vel,mass));
	}
}

void myDisplay(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0,0,400,400);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(-1,1,-1,1,1,-1);

	glClearColor(0,0,0,0);
	glColor3f(1.0,1.0,1.0);

	update_particles();
	
	Particle *temp_particle;
	for (int i = 0; i<NUM_PARTICLES; i++){
		temp_particle = PARTICLES[i];
		glBegin(GL_POINTS);
		glVertex3f(temp_particle->position.x,temp_particle->position.y,temp_particle->position.z);
		glEnd();
	}

	glFlush();
	glutSwapBuffers();
}

int main(int argc, char* argv[]){
	//assuming kernel type gaussian.

	/*
	Simulation involves:
	1. compute density at each particles.
	2. compute pressure from density.
	3. compute forces from pressure.
	4. compute external forces.
	5. Apply those forces to move particles.
	6. Update particle position. (draw as spheres I think in OpenGL)

	Step 1:
	rho_s(r_j) = sum(kernel(r_ij)) (summing over all particles.

	interpreting as r_ij and distance from r_j to r_i. Could be vector as well.

	Step 2:
	pseudo pressure at particle i is:
	p_i = k*(rho_i - rho_0)

	rho_0 is the target density, k is the stiffness constant.

	k = c_s^2 where c_s is speed of sound (in water for us).

	alternate form on intel notes for low compressible fluids. Could use canonical form for simplicity now.

	Step 3:
	pressure gradient at particle j = rho_j*sum(p_i/rho_i^2 + p_j/rho_j^2)*deriv_of_kernel(r_ij)

	Step 4:
	*/

	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

	//The size and position of the window
	glutInitWindowSize(400,400);
	glutInitWindowPosition(0,0);
	glutCreateWindow("Tyler and Zack AS3");

	initScene();

	glutDisplayFunc(myDisplay);
	glutIdleFunc(myDisplay);


	glutMainLoop();

	return 1;
}