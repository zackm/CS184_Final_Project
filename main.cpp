#include "Particle.h"

#include "Triangle.h"

#include <vector>

#include <GL/glut.h>
#include <GL/glu.h>

#include <stdlib.h>
#include <iostream>

#include <cmath>

#include "Vec3.h"

#include "Container.h"

#include "Neighbor.h"

using namespace std;

/*******************
* GLOBAL VARIABLES *
*******************/

Container CONTAINER(Vec3(1,1,1),Vec3(0,0,0));//very simple cube for now. Later, make it a particle itself.
vector<Particle*> PARTICLES;//particles that we do SPH on.
vector<Triangle*> TRIANGLES;//triangles made from marching cubes to render
vector<vector<vector<float> > > GRID_DENSITY;//Grid for marching squares. Probably a better data structure we can use.
vector<vector<vector<bool> > > GRID_BOOL; //bools corresponding to that grid
vector<vector<vector<Vec3> > > VERTEX_MATRIX;//list of vertices corresponding to the densities on the grid.

const float TIMESTEP = .01;//time elapsed between iterations
int NUM_PARTICLES = 350;
Vec3 GRAVITY(0,-9.8f,0);
const float MAX_DISPLACEMENT = 10.0f; //the maximum amount that a particle can be displaced, use when we get NaN.
const float MASS = .01f;//could set it to any number really.
const float STIFFNESS = 0.0f;//10000.0f;//1000.0f;//BOLTZMANN*TEMPERATURE/MOLAR_MASS;// for water; //0 for liquids actually, assuming incompressible.
const float VISCOSITY = 1.004f;// for water;
const float ELECTRIC_CONST = 1e5;
const float DRAG = 1200.0f;
const float SURFACE_TENSION = 1e2;
const float PRESSURE = 10.0f;

const float CUBE_TOL = .125f;//either grid size or tolerance for adaptive cubes, reciprocal must be an integer for now.
const float DENSITY_TOL = 1.5f;//also used for marching grid, for density of the particles

Neighbor NEIGHBOR; //neighbor object used for calculations
const float SUPPORT_RADIUS = .25f;//radius of support used by neighbor function to divide space into grid
const float MAX_KERNEL_RADIUS = SUPPORT_RADIUS/2.0f;

bool USE_ADAPTIVE = false; //for adaptive or uniform marching cubes.

const float PI = 3.1415926;

/*
simple dot product between two vectors.
*/
float dot(Vec3 v1, Vec3 v2){
	return v1.x*v2.x+v1.y*v2.y+v1.z+v2.z;
}

/******************
* Kernel Function *
******************/
//default is used for all terms except pressure and viscosity
float default_kernel(Vec3 r_i, Vec3 r_j){

	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float h = MAX_KERNEL_RADIUS;

	if (sqrt(mag)>h){
		return 0;
	}

	return (315/(64*PI*pow(h,9.0f)))*pow((h*h - mag),3.0f);
}

Vec3 default_gradient(Vec3 r_i,Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float h = MAX_KERNEL_RADIUS;

	Vec3 v(0,0,0);
	if (sqrt(mag)>h){
		return v;
	}

	float coeff = -945.0f/(32.0f*PI*pow(h,9.0f));

	return diff_vec * coeff * pow((h*h - mag),2.0f);
}

float default_laplacian(Vec3 r_i,Vec3 r_j){

	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float h = MAX_KERNEL_RADIUS;

	if (sqrt(mag)>h){
		return 0;
	}

	float coeff = -945.0f/(32.0f*PI*pow(h,9.0f));

	return coeff * (h*h - mag) * (3.0f*h*h - 7.0f*mag);
}

Vec3 pressure_kernel_gradient(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float h = MAX_KERNEL_RADIUS;

	Vec3 v(0,0,0);
	if (sqrt(mag)>h){
		return v;
	}

	float coeff = (45/(PI*pow(h,6.0f)))*pow((h-sqrt(mag)),2.0f);

	return diff_vec*(-coeff)/(sqrt(mag));
}

float viscosity_kernel_laplacian(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float h = MAX_KERNEL_RADIUS;

	if (sqrt(mag)>h){
		return 0;
	}

	return (45.0f/(PI*pow(h,6)))*(h-sqrt(mag));
}

/********************
* Physics Functions *
********************/
Vec3 kinematic_polynomial(Vec3 pos, Vec3 vel, Vec3 acc,float t){
	//this method for running the time step may not be numerically stable enough.
	return (acc*t*t*.5f)+(vel*t)+pos;
}

/*******************
* Particle Methods *
*******************/
float density_at_particle(Particle* part){
	float density = default_kernel(part->position,part->position);
	Particle* temp_particle;
	vector<int> neighbor_vec = part->neighbors;
	for (int i = 0; i<neighbor_vec.size(); i++){ // changed to neighbors

		temp_particle = PARTICLES[neighbor_vec[i]];
		if (i == neighbor_vec[i]) {
			continue;
		}

		density += temp_particle->mass*default_kernel(part->position,temp_particle->position);
	}
	return density;
}

float density_at_point(Vec3 point){
	//first should generate a list of the particles we need to check, using the box for this point.
	//int box_number = NEIGHBOR.compute_box_num(point,SUPPORT_RADIUS,CONTAINER.min.x,CONTAINER.max.x);

	Particle* temp_particle;
	//vector<int> neighbor_vec = NEIGHBOR.box_particles[box_number];
	float density = default_kernel(point,point);
	for (int i = 0; i<PARTICLES.size(); i++){
		//update density. There will be a problem if the point is exactly equal to sum particle (divide by zero error).
		temp_particle = PARTICLES[i];
		//if (i == neighbor_vec[i]) {
		//	continue;
		//}
		density += temp_particle->mass*default_kernel(point,temp_particle->position);
	}

	return density;
}

/*
Create a new list of particles from the old list. Then, throw the old list.
To do this, calculate all quanities in Navier-Stokes, then use timestep to
update particle location from old location and velocity.
*/
void run_time_step(){
	vector<Particle*> new_particles;
	vector<float> density_list;
	vector<float> pressure_list;
	vector<Vec3> pressure_grad_list;
	vector<Vec3> viscosity_list;
	vector<float> color_list;
	//vector<Vec3> normal_list;
	vector<Vec3> tension_list;
	vector<Vec3> repulsion_list;

	//update using slow algorithm for now
	Particle *base_particle, *temp_particle, *new_particle;
	float number_density = 0, density = 0, pressure = 0;

	//Sets density at each point
	for (int i = 0; i<NUM_PARTICLES; i++){
		density = 0;
		base_particle = PARTICLES[i];
		density = density_at_particle(base_particle);
		density_list.push_back(density);

		//changed this to not include stiffness at all
		pressure_list.push_back(STIFFNESS*(PRESSURE));
	}

	//Sets pressure gradient at each point using densities from last loop
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 pressure_gradient(0,0,0);

		vector<int> neighbor_vec = base_particle->neighbors;
		int n = neighbor_vec.size();
		for (int j = 0; j<n; j++){ // changed to neighbors
			temp_particle = PARTICLES[neighbor_vec[j]];
			if (i == neighbor_vec[j]) {
				continue;
			}

			Vec3 weight = pressure_kernel_gradient(base_particle->position,temp_particle->position);
			pressure_gradient += weight * temp_particle->mass * ((pressure_list[i]+pressure_list[j])/(2.0f*density_list[j])); 

		}
		pressure_grad_list.push_back(pressure_gradient);
	}

	//Sets viscosity at each particle
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 viscosity_laplacian(0,0,0);
		vector<int> neighbor_vec = base_particle->neighbors;
		for (int j = 0; j<neighbor_vec.size(); j++){ // changed to neighbors
			temp_particle = PARTICLES[neighbor_vec[j]];

			float weight = viscosity_kernel_laplacian(base_particle->position,temp_particle->position);
			viscosity_laplacian += ((temp_particle->velocity - base_particle->velocity)/density_list[i])*weight * temp_particle->mass;
		}
		viscosity_list.push_back(viscosity_laplacian*VISCOSITY);
	}

	//Sets color field at each point for surface tension
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		float color = 0.0f;
		vector<int> neighbor_vec = base_particle->neighbors;
		for (int j = 0; j<neighbor_vec.size(); j++){
			temp_particle = PARTICLES[neighbor_vec[j]];
			color += (temp_particle->mass / density_list[j]) * default_laplacian(base_particle->position,temp_particle->position);

		}
		color_list.push_back(color);
	}

	//set the normal at each point
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 normal(0,0,0);
		vector<int> neighbor_vec = base_particle->neighbors;
		for (int j = 0; j<neighbor_vec.size(); j++){
			temp_particle = PARTICLES[neighbor_vec[j]];
			normal += default_gradient(base_particle->position,temp_particle->position)*(temp_particle->mass / density_list[j]);

		}

		float length = sqrt(dot(normal,normal));
		if(length>0.0){
			normal = normal/sqrt(length);
			tension_list.push_back(normal*(-1.0f*color_list[i]));
		}else{
			Vec3 v(0,0,0);
			tension_list.push_back(v);
		}
	}

	//set the repulsion force at each point
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 repulsion(0,0,0);
		vector<int> neighbor_vec = base_particle->neighbors;
		for (int j = 0; j<neighbor_vec.size(); j++){
			temp_particle = PARTICLES[neighbor_vec[j]];
			Vec3 diff_vec = temp_particle->position - base_particle->position;
			float length = sqrt(dot(diff_vec,diff_vec));
			Vec3 electric_repulsion;
			electric_repulsion = diff_vec/(max(pow(length,3.0f),.01f));
			
			repulsion += electric_repulsion * (default_kernel(base_particle->position,temp_particle->position))*(temp_particle->mass / density_list[j]);

		}
		repulsion_list.push_back(repulsion);
	}

	//Create new particles from old particles and from pressure gradient.
	for (int i = 0; i<NUM_PARTICLES; i++){
		temp_particle = PARTICLES[i];

		//using Navier Stokes, calculate the change in velocity.
		Vec3 velocity = temp_particle->velocity;

		//First add external forces
		Vec3 acceleration = GRAVITY - ((velocity*DRAG) + tension_list[i]*SURFACE_TENSION
						    + (viscosity_list[i]*VISCOSITY)
							+ repulsion_list[i]*ELECTRIC_CONST
							- pressure_grad_list[i])/density_list[i];

		Vec3 position = temp_particle->position;

		Vec3 new_position = kinematic_polynomial(position,velocity,acceleration,TIMESTEP);
		Vec3 new_velocity = temp_particle->velocity + acceleration*TIMESTEP;
		float mass = temp_particle->mass;

		temp_particle = new Particle(new_position,new_velocity,mass);

		CONTAINER.in_container(temp_particle,TIMESTEP); //applies reflections if outside of boundary.

		new_particles.push_back(temp_particle);
	}

	PARTICLES = new_particles;

	//reset neighbor structure 
	NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER);
}



/*****************
* OpenGL Methods *
*****************/
void initScene(){
	//create a list of random particles
	Particle* new_part;
	float noise = .5f;
	float x,y,z,v_x,v_y,v_z;

	////Semi random grid of particles
	//float step = .05f;
	//for(float i = CONTAINER.min.x; i<(CONTAINER.max.x); i=i+step){
	//	for(float j = CONTAINER.min.y; j<(CONTAINER.max.y/2.0f); j=j+step){
	//		noise = float(rand())/(float(RAND_MAX))*.05f;
	//		Vec3 pos(i+noise,j+noise,0);
	//		Vec3 vel(0,0,0);
	//		new_part = new Particle(pos,vel,MASS);
	//		PARTICLES.push_back(new_part);
	//	}
	//}
	//NUM_PARTICLES = PARTICLES.size();

	//random particles
	for (int i = 0; i<NUM_PARTICLES; i++){
		x = .2f+float(rand())/(float(RAND_MAX))*.1f;
		y = float(rand())/(float(RAND_MAX))/5.0 + .4f;
		z = 0.0f;//.2f + float(rand())/(float(RAND_MAX))*.1f;

		v_x = -.5f+float(rand())/(float(RAND_MAX))*2.0f*noise;
		v_y = -0.2+float(rand())/(float(RAND_MAX))*noise;
		v_z = 0.0f;//-0.2+float(rand())/(float(RAND_MAX))*noise;


		Vec3 pos(x,y,z);
		Vec3 vel(v_x,v_y,v_z);
		float mass = 4.0f+(float(rand())/(float(RAND_MAX)))*5.0f;
		//also need to instantiate the other fields
		new_part = new Particle(pos,vel,MASS);
		PARTICLES.push_back(new_part);
	}
	NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER);

	glEnable(GL_DEPTH_TEST);

}

void myDisplay(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0,0,400,400);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(90,1.0f,.1,-1000);
	glOrtho(CONTAINER.min.x,CONTAINER.max.x,CONTAINER.min.y,CONTAINER.max.y,CONTAINER.min.z,CONTAINER.max.z);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//gluLookAt(1.2,1.2,1.2,0,0,0,0,1,0);

	run_time_step();

	//marching_cubes();

	//draw particles
	Particle* temp_part;
	for (int i = 0; i<PARTICLES.size(); i++){
		temp_part = PARTICLES[i];
		glClearColor(0,0,0,0);
		glPointSize(4.0f);
		glBegin(GL_POINTS);
		// alternate particle colors depending on box in grid
		/* if (temp_part->box % 2 == 0) {
		glColor3f(1.0,0,0);
		} else {
		glColor3f(0,1.0,1.0);
		}*/
		glColor3f(0,0,1.0);
		glVertex3f(temp_part->position.x,temp_part->position.y,temp_part->position.z);
		glEnd();
	}

	////draw triangles
	//Triangle *temp_triangle;
	//for (int i = 0; i<TRIANGLES.size(); i++){
	//	temp_triangle = TRIANGLES[i];

	//	glClearColor(0,0,0,0);
	//	glColor3f(0,0,1.0f);
	//	glBegin(GL_TRIANGLES);
	//	glVertex3f(temp_triangle->a.x,temp_triangle->a.y,temp_triangle->a.z);
	//	glVertex3f(temp_triangle->b.x,temp_triangle->b.y,temp_triangle->b.z);
	//	glVertex3f(temp_triangle->c.x,temp_triangle->c.y,temp_triangle->c.z);
	//	glEnd();
	//}

	////Draw wireframe container
	//glPolygonMode(GL_FRONT, GL_LINE);
	//glPolygonMode(GL_BACK, GL_LINE);

	//glDisable(GL_LIGHTING);
	//glClearColor (0.0, 0.0, 0.0, 0.0);
	//glColor3f(1.0f,1.0f,1.0f);

	//glBegin(GL_POLYGON);
	//glVertex3f(CONTAINER.min.x,CONTAINER.min.y,CONTAINER.max.z);
	//glVertex3f(CONTAINER.min.x,CONTAINER.max.y,CONTAINER.max.z);
	//glVertex3f(CONTAINER.max.x,CONTAINER.max.y,CONTAINER.max.z);
	//glVertex3f(CONTAINER.max.x,CONTAINER.min.y,CONTAINER.max.z);
	//glEnd();

	//glBegin(GL_POLYGON);
	//glVertex3f(CONTAINER.min.x,CONTAINER.min.y,CONTAINER.min.z);
	//glVertex3f(CONTAINER.min.x,CONTAINER.max.y,CONTAINER.min.z);
	//glVertex3f(CONTAINER.max.x,CONTAINER.max.y,CONTAINER.min.z);
	//glVertex3f(CONTAINER.max.x,CONTAINER.min.y,CONTAINER.min.z);
	//glEnd();

	//glBegin(GL_LINES);
	//glVertex3f(CONTAINER.min.x,CONTAINER.min.y,CONTAINER.min.z);
	//glVertex3f(CONTAINER.min.x,CONTAINER.min.y,CONTAINER.max.z);
	//glEnd();

	//glBegin(GL_LINES);
	//glVertex3f(CONTAINER.min.x,CONTAINER.max.y,CONTAINER.min.z);
	//glVertex3f(CONTAINER.min.x,CONTAINER.max.y,CONTAINER.max.z);
	//glEnd();

	//glBegin(GL_LINES);
	//glVertex3f(CONTAINER.max.x,CONTAINER.max.y,CONTAINER.min.z);
	//glVertex3f(CONTAINER.max.x,CONTAINER.max.y,CONTAINER.max.z);
	//glEnd();

	//glBegin(GL_LINES);
	//glVertex3f(CONTAINER.max.x,CONTAINER.min.y,CONTAINER.min.z);
	//glVertex3f(CONTAINER.max.x,CONTAINER.min.y,CONTAINER.max.z);
	//glEnd();

	//glPolygonMode(GL_FRONT, GL_FILL); // fill mode
	//glPolygonMode(GL_BACK, GL_FILL);

	glPopMatrix();

	glFlush();
	glutSwapBuffers();
}

int main(int argc, char* argv[]){
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

	//The size and position of the window
	glutInitWindowSize(400,400);
	glutInitWindowPosition(0,0);
	glutCreateWindow("Tyler and Zack Final Project");

	initScene();

	glutDisplayFunc(myDisplay);
	glutIdleFunc(myDisplay);

	glutMainLoop();

	return 0;
}