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
const float LIFETIME = 100.0f;
float CURRENT_TIME = 0.0f;
int NUM_PARTICLES = 0;
Vec3 GRAVITY(0,-9.8f,0);
const float MASS = .02f;//could set it to any number really.
const float IDEAL_DENSITY = 998.0f;
const float STIFFNESS = 3.0f;//for pressure difference
const float VISCOSITY = 3.5f;
const float SURFACE_TENSION = .0728f;
const float KERNEL_PARTICLES = 20;
const float TENSION_THRESHOLD = sqrt(IDEAL_DENSITY/KERNEL_PARTICLES);

const float CUBE_TOL = .125f;//either grid size or tolerance for adaptive cubes, reciprocal must be an integer for now.
const float DENSITY_TOL = 1.5f;//also used for marching grid, for density of the particles

Neighbor NEIGHBOR; //neighbor object used for calculations
const float H = .26;
const float SUPPORT_RADIUS = .05;//2*H;


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

	return (315/(64*PI*pow(H,9.0f)))*pow((H*H - mag),3.0f);
}

Vec3 default_gradient(Vec3 r_i,Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float coeff = -945.0f/(32.0f*PI*pow(H,9.0f));

	return diff_vec * coeff * pow((H*H - mag),2.0f);
}

float default_laplacian(Vec3 r_i,Vec3 r_j){

	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float coeff = -945.0f/(32.0f*PI*pow(H,9.0f));

	return coeff * (H*H - mag) * (3.0f*H*H - 7.0f*mag);
}

Vec3 pressure_kernel_gradient(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float coeff = (45/(PI*pow(H,6.0f)))*pow((H-sqrt(mag)),2.0f);

	Vec3 v(0,0,0);
	return diff_vec*(-coeff)/(sqrt(mag));
}

float viscosity_kernel_laplacian(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	return (45.0f/(PI*pow(H,6)))*(H-sqrt(mag));
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
float density_at_point(Vec3 point){
	//first should generate a list of the particles we need to check, using the box for this point.
	//int box_number = NEIGHBOR.compute_box_num(point,SUPPORT_RADIUS,CONTAINER.min.x,CONTAINER.max.x);

	Particle* temp_particle;
	//vector<int> neighbor_vec = NEIGHBOR.box_particles[box_number];
	float density = default_kernel(point,point);
	for (int i = 0; i<PARTICLES.size(); i++){
		//update density. There will be a problem if the point is exactly equal to sum particle (divide by zero error).
		temp_particle = PARTICLES[i];

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
	vector<Vec3> tension_list;

	//update using slow algorithm for now
	Particle *base_particle, *temp_particle, *new_particle;
	float density = 0;

	//Sets density at each point
	for (int i = 0; i<NUM_PARTICLES; i++){
		density = 0;
		base_particle = PARTICLES[i];


		for (int j = 0; j<NUM_PARTICLES; j++){ // changed to neighbors
			if(i!=j){
				temp_particle = PARTICLES[i];
				Vec3 r = base_particle->position-temp_particle->position;
				float mag = dot(r,r);
				if(mag<H*H){
					density += temp_particle->mass*default_kernel(base_particle->position,temp_particle->position);
				}
			}
		}
		density_list.push_back(density);

		//changed this to not include stiffness at all
		pressure_list.push_back(STIFFNESS*(density-IDEAL_DENSITY));
	}

	//Sets pressure gradient at each point using densities from last loop
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 pressure_gradient(0,0,0);

		vector<int> neighbor_vec = base_particle->neighbors;
		int n = neighbor_vec.size();
		for (int j = 0; j<NUM_PARTICLES; j++){ // changed to neighbors
			if(i!=j){
				temp_particle = PARTICLES[j];
				Vec3 r = base_particle->position-temp_particle->position;
				float mag = dot(r,r);

				if(mag<H*H){
					Vec3 weight = pressure_kernel_gradient(base_particle->position,temp_particle->position);
					pressure_gradient += weight * temp_particle->mass * ((pressure_list[i]+pressure_list[j])/(2.0f*density_list[j])); 
				}
			}

			//if(mag==0){
			//	cout<<'h'<<endl;
			//}
		}
		pressure_grad_list.push_back(pressure_gradient*(-1.0f));
	}

	//Sets viscosity at each particle
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 viscosity_laplacian(0,0,0);
		vector<int> neighbor_vec = base_particle->neighbors;
		for (int j = 0; j<NUM_PARTICLES; j++){ // changed to neighbors
			if(i!=j){
				temp_particle = PARTICLES[j];
				Vec3 r = base_particle->position-temp_particle->position;
				float mag = dot(r,r);
				if(mag<H*H){
					float weight = viscosity_kernel_laplacian(base_particle->position,temp_particle->position);
					viscosity_laplacian += ((temp_particle->velocity - base_particle->velocity)/density_list[i])*weight * temp_particle->mass;
				}
			}
		}
		viscosity_list.push_back(viscosity_laplacian*VISCOSITY);
	}

	//Sets color field laplacian at each point for surface tension
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		float color = 0.0f;
		vector<int> neighbor_vec = base_particle->neighbors;
		for (int j = 0; j<NUM_PARTICLES; j++){
			if(i!=j){
				temp_particle = PARTICLES[j];
				Vec3 r = base_particle->position-temp_particle->position;
				float mag = dot(r,r);
				if(mag<H*H){
					color += (temp_particle->mass / density_list[j]) * default_laplacian(base_particle->position,temp_particle->position);
				}
			}
		}
		color_list.push_back(color);
	}

	//set the normal at each point
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 normal(0,0,0);
		vector<int> neighbor_vec = base_particle->neighbors;
		for (int j = 0; j<NUM_PARTICLES; j++){
			if(i!=j){
				temp_particle = PARTICLES[j];
				Vec3 r = base_particle->position-temp_particle->position;
				float mag = dot(r,r);
				if(mag<H*H){
					normal += default_gradient(base_particle->position,temp_particle->position)*(temp_particle->mass / density_list[j]);
				}
			}

		}

		float length = sqrt(dot(normal,normal));
		//cout<<length<<endl;
		if(length>TENSION_THRESHOLD){
			//cout<<'n'<<endl;
			normal = normal/sqrt(length);
			tension_list.push_back(normal*(-1.0f*color_list[i]));
		}else{
			Vec3 v(0,0,0);
			tension_list.push_back(v);
		}
	}


	//Create new particles from old particles and from pressure gradient.
	for (int i = 0; i<NUM_PARTICLES; i++){
		temp_particle = PARTICLES[i];

		//using Navier Stokes, calculate the change in velocity.
		Vec3 velocity = temp_particle->velocity;

		//First add external forces
		Vec3 acceleration = GRAVITY + (tension_list[i]*SURFACE_TENSION
			+ (viscosity_list[i]*VISCOSITY)
			+ pressure_grad_list[i])/density_list[i];

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
	//NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER);
}

/*****************
* OpenGL Methods *
*****************/
void initScene(){
	//create a list of random particles
	Particle* new_part;
	float noise = float(rand())/(float(RAND_MAX))*.1;
	float x,y,z,v_x,v_y,v_z;

	////Semi random grid of particles
	//float step = .025;
	//for(float i = 2.0*CONTAINER.max.x/5.0f; i<3.0f*(CONTAINER.max.x)/5.0f; i=i+step){
	//	for(float j = 3.0*CONTAINER.max.y/4.0f; j<(CONTAINER.max.y); j=j+step){
	//		noise = float(rand())/(float(RAND_MAX))*.05f;
	//		Vec3 pos(i+float(rand())/(float(RAND_MAX))*.1,j+float(rand())/(float(RAND_MAX))*.1,0);
	//		Vec3 vel(0,0,0);
	//		new_part = new Particle(pos,vel,MASS);
	//		PARTICLES.push_back(new_part);
	//	}
	//}

	//float step = .005;
	//for(float i = CONTAINER.min.x; i<(CONTAINER.max.x); i=i+step){
	//	for(float j = CONTAINER.min.y; j<(CONTAINER.max.y)/2.0f; j=j+step){
	//		noise = float(rand())/(float(RAND_MAX))*.05f;
	//		Vec3 pos(i,j,0);
	//		Vec3 vel(0,0,0);
	//		new_part = new Particle(pos,vel,MASS);
	//		PARTICLES.push_back(new_part);
	//	}
	//}
	//NUM_PARTICLES = PARTICLES.size();

	////random particles
	//for (int i = 0; i<NUM_PARTICLES; i++){
	//	x = .2f+float(rand())/(float(RAND_MAX))*.1f;
	//	y = float(rand())/(float(RAND_MAX))/5.0 + .1f;
	//	z = 0.0f;//.2f + float(rand())/(float(RAND_MAX))*.1f;

	//	v_x = -.5f+float(rand())/(float(RAND_MAX))*2.0f*noise;
	//	v_y = -0.2+float(rand())/(float(RAND_MAX))*noise;
	//	v_z = 0.0f;//-0.2+float(rand())/(float(RAND_MAX))*noise;


	//	Vec3 pos(x,y,z);
	//	Vec3 vel(v_x,v_y,v_z);
	//	float mass = 4.0f+(float(rand())/(float(RAND_MAX)))*5.0f;
	//	//also need to instantiate the other fields
	//	new_part = new Particle(pos,vel,MASS);
	//	PARTICLES.push_back(new_part);
	//}
	////NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER);

	//glEnable(GL_DEPTH_TEST);

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
	CURRENT_TIME += TIMESTEP;
	//marching_cubes();

	if(CURRENT_TIME>LIFETIME){
		exit(0);
	}

	if(CURRENT_TIME<2.2){
		//throw in a new particle.
		float noise = float(rand())/(float(RAND_MAX))*.05f;
		Vec3 pos(.1+noise,.9,0);
		Vec3 vel(float(rand())/(float(RAND_MAX)),float(rand())/(float(RAND_MAX)),0);
		Particle* new_part = new Particle(pos,vel,MASS);
		PARTICLES.push_back(new_part);
		NUM_PARTICLES++;
	}
	NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER);
	
	//draw particles
	glPointSize(4.0f);
	glBegin(GL_POINTS);
	Particle* temp_part;
	for (int i = 0; i<PARTICLES.size(); i++){
		temp_part = PARTICLES[i];
		glClearColor(0,0,0,0);


		// alternate particle colors depending on box in grid
		/*if (temp_part->box % 2 == 0) {
		glColor3f(1.0,0,0);
		} else {
		glColor3f(0,1.0,1.0);
		}*/
		glColor3f(0,0,1.0);
		glVertex3f(temp_part->position.x,temp_part->position.y,temp_part->position.z);
        
//        glColor3f(1,0,0);
//        glPushMatrix();
//        glTranslated(temp_part->position.x,temp_part->position.y,temp_part->position.z);
//        glutWireSphere(H,10,10);
//        glPopMatrix();
	}
	glEnd();

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

	//Draw wireframe container
	glPolygonMode(GL_FRONT, GL_LINE);
	glPolygonMode(GL_BACK, GL_LINE);

	glDisable(GL_LIGHTING);
	glClearColor (0.0, 0.0, 0.0, 0.0);
	glColor3f(1.0f,1.0f,1.0f);

	glBegin(GL_POLYGON);
	glVertex3f(CONTAINER.min.x,CONTAINER.min.y,CONTAINER.max.z);
	glVertex3f(CONTAINER.min.x,CONTAINER.max.y,CONTAINER.max.z);
	glVertex3f(CONTAINER.max.x,CONTAINER.max.y,CONTAINER.max.z);
	glVertex3f(CONTAINER.max.x,CONTAINER.min.y,CONTAINER.max.z);
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3f(CONTAINER.min.x,CONTAINER.min.y,CONTAINER.min.z);
	glVertex3f(CONTAINER.min.x,CONTAINER.max.y,CONTAINER.min.z);
	glVertex3f(CONTAINER.max.x,CONTAINER.max.y,CONTAINER.min.z);
	glVertex3f(CONTAINER.max.x,CONTAINER.min.y,CONTAINER.min.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(CONTAINER.min.x,CONTAINER.min.y,CONTAINER.min.z);
	glVertex3f(CONTAINER.min.x,CONTAINER.min.y,CONTAINER.max.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(CONTAINER.min.x,CONTAINER.max.y,CONTAINER.min.z);
	glVertex3f(CONTAINER.min.x,CONTAINER.max.y,CONTAINER.max.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(CONTAINER.max.x,CONTAINER.max.y,CONTAINER.min.z);
	glVertex3f(CONTAINER.max.x,CONTAINER.max.y,CONTAINER.max.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(CONTAINER.max.x,CONTAINER.min.y,CONTAINER.min.z);
	glVertex3f(CONTAINER.max.x,CONTAINER.min.y,CONTAINER.max.z);
	glEnd();

	glPolygonMode(GL_FRONT, GL_FILL); // fill mode
	glPolygonMode(GL_BACK, GL_FILL);

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



//float density_at_particle(Particle* part){
//	float density = default_kernel(part->position,part->position);
//	Particle* temp_particle;
//	vector<int> neighbor_vec = part->neighbors;
//
//	for (int i = 0; i<NUM_PARTICLES; i++){ // changed to neighbors
//		temp_particle = PARTICLES[i];
//		density += temp_particle->mass*default_kernel(part->position,temp_particle->position);
//	}
//	//cout<<density<<endl;
//	return density;
//}
