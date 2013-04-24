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
Container CONTAINER(Vec3(.5,.5,.5),Vec3(0,0,0));//very simple cube for now. Later, make it a particle itself.
vector<Particle*> PARTICLES;//particles that we do SPH on.
vector<Triangle*> TRIANGLES;//triangles made from marching cubes to render
vector<vector<vector<float> > > GRID_DENSITY;//Grid for marching squares. Probably a better data structure we can use.
vector<vector<vector<bool> > > GRID_BOOL; //bools corresponding to that grid
vector<vector<vector<Vec3> > > VERTEX_MATRIX;//list of vertices corresponding to the densities on the grid.

const float TIMESTEP = .005;//time elapsed between iterations
const float LIFETIME = 100.0f;
float CURRENT_TIME = 0.0f;
int NUM_PARTICLES = 0;
Vec3 GRAVITY(0,-9.8,0);
const float MASS = .02f;//could set it to any number really.
const float IDEAL_DENSITY = 1000.0f;
const float STIFFNESS = 3.0f;//for pressure difference
const float VISCOSITY = 3.5f;
const float SURFACE_TENSION = .07f;
const float TENSION_THRESHOLD = 7.0f;

const float CUBE_TOL = .02f;//either grid size or tolerance for adaptive cubes, reciprocal must be an integer for now.
const float DENSITY_TOL = 1000.0f;//also used for marching grid, for density of the particles

Neighbor NEIGHBOR; //neighbor object used for calculations
const float H = .05;
const float SUPPORT_RADIUS = .1;

bool RENDERING_ISOSURFACE = false;
bool USE_ADAPTIVE = false; //for adaptive or uniform marching cubes.

const float PI = 3.1415926;
const float DRAW_RADIUS = .01f;

void marching_cubes();

/*
simple dot product between two vectors.
*/
float dot(Vec3 v1, Vec3 v2){
	return (v1.x*v2.x)+(v1.y*v2.y)+(v1.z*v2.z);
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
	int box_number = NEIGHBOR.compute_box_num(point,SUPPORT_RADIUS,CONTAINER.min.x,CONTAINER.max.x);

	Particle* temp_particle;
	vector<int> neighbor_vec = NEIGHBOR.box_particles[box_number];
	float density = 0;//default_kernel(point,point);
	for (int i = 0; i<neighbor_vec.size(); i++){
		//update density. There will be a problem if the point is exactly equal to sum particle (divide by zero error).
		temp_particle = PARTICLES[neighbor_vec[i]];
		Vec3 diff_vec = point-temp_particle->position;
		float mag = dot(diff_vec,diff_vec);
		if(mag<H*H){
			density += temp_particle->mass*default_kernel(point,temp_particle->position);
		}
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

	NEIGHBOR.place_particles(PARTICLES, SUPPORT_RADIUS, CONTAINER);

	//update using slow algorithm for now
	Particle *base_particle, *temp_particle, *new_particle;
	float density = 0;

	//Sets density at each point
	for (int i = 0; i<NUM_PARTICLES; i++){
		density = 0;
		base_particle = PARTICLES[i];

		vector<int> neighbor_vec = base_particle->neighbors;
		int n = neighbor_vec.size();
		for (int j = 0; j<n; j++){ // changed to neighbors
			temp_particle = PARTICLES[neighbor_vec[j]];
			Vec3 r = base_particle->position-temp_particle->position;
			float mag = dot(r,r);
			if(mag<H*H){
				density += temp_particle->mass*default_kernel(base_particle->position,temp_particle->position);
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
		for (int j = 0; j<neighbor_vec.size(); j++){ // changed to neighbors
			if(i!=neighbor_vec[j]){
				temp_particle = PARTICLES[neighbor_vec[j]];
				Vec3 r = base_particle->position-temp_particle->position;
				float mag = dot(r,r);

				if(mag<H*H){
					Vec3 weight = pressure_kernel_gradient(base_particle->position,temp_particle->position);
					pressure_gradient += weight * temp_particle->mass * ((pressure_list[i]+pressure_list[j])/(2.0f*density_list[j])); 
				}
			}
		}
		pressure_grad_list.push_back(pressure_gradient*(-1.0f));
	}

	//Sets viscosity at each particle
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 viscosity_laplacian(0,0,0);
		vector<int> neighbor_vec = base_particle->neighbors;
		for (int j = 0; j<neighbor_vec.size(); j++){
			if(i!=neighbor_vec[j]){
				temp_particle = PARTICLES[neighbor_vec[j]];
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
		for (int j = 0; j<neighbor_vec.size(); j++){
			if(i!=neighbor_vec[j]){
				temp_particle = PARTICLES[neighbor_vec[j]];
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
		for (int j = 0; j<neighbor_vec.size(); j++){
			if(i!=neighbor_vec[j]){
				temp_particle = PARTICLES[neighbor_vec[j]];
				Vec3 r = base_particle->position-temp_particle->position;
				float mag = dot(r,r);
				if(mag<H*H){
					normal += default_gradient(base_particle->position,temp_particle->position)*(temp_particle->mass / density_list[j]);
				}
			}

		}

		float length = sqrt(dot(normal,normal));
		if(length>TENSION_THRESHOLD){
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
	NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER);
}

/*****************
* OpenGL Methods *
*****************/
void initScene(){
	//create a list of random particles
	Particle* new_part;
	float noise = float(rand())/(float(RAND_MAX))*.1;
	float x,y,z,v_x,v_y,v_z;

	//2D Drop Scene
	float step = .017;
	for(float i = 2.0*CONTAINER.max.x/5.0f; i<3.0*(CONTAINER.max.x)/5.0f; i=i+step){
		for(float j = 4.0*CONTAINER.max.y/5.0f; j<(CONTAINER.max.y); j=j+step){
			noise = float(rand())/(float(RAND_MAX))*.05f;
			Vec3 pos(i,j,0);
			Vec3 vel(0,-5,0);
			new_part = new Particle(pos,vel,MASS);
			PARTICLES.push_back(new_part);
		}
	}

	step = .017;
	for(float i = CONTAINER.min.x; i<(CONTAINER.max.x); i=i+step){
		for(float j = CONTAINER.min.y; j<(CONTAINER.max.y)/5.0f; j=j+step){
			noise = float(rand())/(float(RAND_MAX))*.05f;
			Vec3 pos(i,j,0);
			Vec3 vel(0,0,0);
			new_part = new Particle(pos,vel,MASS);
			PARTICLES.push_back(new_part);
		}
	}

	//2D Throw Scene
	//Semi random grid of particles
	/*float step = .007;
	for(float i = 4.0*CONTAINER.max.x/5.0f; i<(CONTAINER.max.x); i=i+step){
	for(float j = 3.0*CONTAINER.max.y/4.0f; j<(CONTAINER.max.y); j=j+step){
	noise = float(rand())/(float(RAND_MAX))*.05f;
	Vec3 pos(i,j,0);
	Vec3 vel(-1,-8,0);
	new_part = new Particle(pos,vel,MASS);
	PARTICLES.push_back(new_part);
	}
	}

	step = .007;
	for(float i = CONTAINER.min.x; i<1.0f*(CONTAINER.max.x)/5.0f; i=i+step){
	for(float j = 3.0*CONTAINER.max.y/4.0f; j<(CONTAINER.max.y); j=j+step){
	noise = float(rand())/(float(RAND_MAX))*.05f;
	Vec3 pos(i,j,0);
	Vec3 vel(5,-5,0);
	new_part = new Particle(pos,vel,MASS);
	PARTICLES.push_back(new_part);
	}
	}*/

	//////3D Drop Scene
	//float step = .02;
	//for(float i = 2.0*CONTAINER.max.x/5.0; i<3.0f*(CONTAINER.max.x)/5.0f; i=i+step){
	//	for(float j = 2.0*CONTAINER.max.y/5.0f; j<3.0f*(CONTAINER.max.y)/5.0f; j=j+step){
	//		for(float k = 1.0*CONTAINER.max.y/5.0f; k<4.0f*(CONTAINER.max.y)/5.0f; k=k+step){
	//			//noise = float(rand())/(float(RAND_MAX))*.05f;
	//			Vec3 pos(i,j,k);
	//			Vec3 vel(0,-3,0);
	//			new_part = new Particle(pos,vel,MASS);
	//			PARTICLES.push_back(new_part);
	//		}
	//	}
	//}

	//step = .02;
	//for(float i = CONTAINER.min.x; i<(CONTAINER.max.x); i=i+step){
	//	for(float j = CONTAINER.min.y; j<1.0f*(CONTAINER.max.y)/5.0f; j=j+step){
	//		for(float k = CONTAINER.min.z; k<(CONTAINER.max.z); k=k+step){
	//			//noise = float(rand())/(float(RAND_MAX))*.05f;
	//			Vec3 pos(i,j,k);
	//			Vec3 vel(0,0,0);
	//			new_part = new Particle(pos,vel,MASS);
	//			PARTICLES.push_back(new_part);
	//		}
	//	}
	//}


	NUM_PARTICLES = PARTICLES.size();
	cout<<NUM_PARTICLES<<endl;
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
	NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER);

	//create some lights
	GLfloat light_position[] = {1,1,1,0};
	GLfloat mat_specular[] = {0,0,0,1.0};
	GLfloat mat_diffuse[] = {1.0,1.0,1.0,1.0};
	GLfloat mat_ambient[] = {.1,.1,.1,1};
	GLfloat mat_shininess[] = {20.0};

	glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
	glMaterialfv(GL_FRONT,GL_SHININESS,mat_shininess);
	glMaterialfv(GL_FRONT,GL_DIFFUSE,mat_diffuse);
	glMaterialfv(GL_FRONT,GL_AMBIENT,mat_ambient);
	glLightfv(GL_LIGHT0,GL_POSITION,light_position);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
}

void myDisplay(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0,0,400,400);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(50,1.0f,.0001,1000);
	//glOrtho(CONTAINER.min.x,CONTAINER.max.x,CONTAINER.min.y,CONTAINER.max.y,CONTAINER.min.z,CONTAINER.max.z);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(.25f,.25f,1.25f,.25f,.25f,0.0f,0,1,0);

	run_time_step();
	CURRENT_TIME += TIMESTEP;

	if(CURRENT_TIME>LIFETIME){
		exit(0);
	}

	//draw particles
	//glPointSize(4.0f);
	//glEnable(GL_LIGHTING);
	if(RENDERING_ISOSURFACE){
		////draw triangles
		marching_cubes();
		Triangle *temp_triangle;
		for (int i = 0; i<TRIANGLES.size(); i++){
			temp_triangle = TRIANGLES[i];

			glClearColor(0,0,0,0);
			glColor3f(0,0,1.0f);
			glBegin(GL_TRIANGLES);
			glVertex3f(temp_triangle->a.x,temp_triangle->a.y,temp_triangle->a.z);
			glVertex3f(temp_triangle->b.x,temp_triangle->b.y,temp_triangle->b.z);
			glVertex3f(temp_triangle->c.x,temp_triangle->c.y,temp_triangle->c.z);
			glEnd();
		}
	}else{
		//draw particles
		glEnable(GL_LIGHTING);
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
			//glColor3f(0,0,1.0);
			//glVertex3f(temp_part->position.x,temp_part->position.y,temp_part->position.z);

			//Draw sphere of radius H around particles
			//glColor3f(0,0,1.0);
			glPushMatrix();
			glTranslated(temp_part->position.x,temp_part->position.y,temp_part->position.z);
			glutSolidSphere(DRAW_RADIUS,16,16);
			glPopMatrix();
		}
	}

	//Draw wireframe container
	glPolygonMode(GL_FRONT, GL_LINE);
	glPolygonMode(GL_BACK, GL_LINE);

	glDisable(GL_LIGHTING);
	glClearColor (0, 0, 0, 0);
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

/************************
* Marching Cube Methods *
************************/
/*
Returns a list of pointers to the different triangles that need to be made for the marching cube case 
corresponding to the input int. Should call some global list of all possible triangles (up to transformation)
to make a particular triangle.
*/
vector<Triangle*> make_triangles(int i, int j, int k){
	//cube define by the corner (i,j), (i,j+1), (i+1,j), and (i+1,j+1).

	int cube_case = GRID_BOOL[i][j][k]*8+GRID_BOOL[i][j+1][k]*4+GRID_BOOL[i+1][j][k]*2+GRID_BOOL[i+1][j+1][k];

	Vec3 vertex_1 = VERTEX_MATRIX[i][j][k];
	Vec3 vertex_2 = VERTEX_MATRIX[i][j+1][k];
	Vec3 vertex_3 = VERTEX_MATRIX[i+1][j][k];
	Vec3 vertex_4 = VERTEX_MATRIX[i+1][j+1][k];
	float weight_1 = GRID_DENSITY[i][j][k];
	float weight_2 = GRID_DENSITY[i][j+1][k];
	float weight_3 = GRID_DENSITY[i+1][j][k];
	float weight_4 = GRID_DENSITY[i+1][j+1][k];

	vector<Triangle*> tri_list;
	Triangle* tri1, *tri2, *tri3;
	Vec3 midpoint_1, midpoint_2;

	if(cube_case==0){
		//no triangles, no corners activated
	}else if(cube_case==1){
		//only corner 11 is turned on.
		tri1 = new Triangle(vertex_4,vertex_2,vertex_3,weight_4,weight_2,weight_3,DENSITY_TOL);
		tri_list.push_back(tri1);
	}else if(cube_case==2){
		//only corner 10 is turned on
		tri1 = new Triangle(vertex_3,vertex_4,vertex_1,weight_3,weight_4,weight_1,DENSITY_TOL);
		tri_list.push_back(tri1);
	}else if(cube_case==3){
		//10 and 11
		midpoint_1 = vertex_3+(vertex_1-vertex_3)*((DENSITY_TOL-weight_3)/(weight_1-weight_3));
		tri1 = new Triangle(vertex_3,vertex_4,midpoint_1);

		midpoint_2 = vertex_4+(vertex_2-vertex_4)*((DENSITY_TOL-weight_4)/(weight_2-weight_4));
		tri2 = new Triangle(vertex_4,midpoint_1,midpoint_2);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
	}else if(cube_case==4){
		//only corner 01 is turned on
		tri1 = new Triangle(vertex_2,vertex_1,vertex_4,weight_2,weight_1,weight_4,DENSITY_TOL);
		tri_list.push_back(tri1);
	}else if(cube_case==5){
		//01 and 11 turned on
		midpoint_1 = vertex_2+(vertex_1-vertex_2)*((DENSITY_TOL-weight_2)/(weight_1-weight_2));
		tri1 = new Triangle(vertex_4,vertex_2,midpoint_1);

		midpoint_2 = vertex_4+(vertex_3-vertex_4)*((DENSITY_TOL-weight_4)/(weight_3-weight_4));
		tri2 = new Triangle(vertex_4,midpoint_1,midpoint_2);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
	}else if(cube_case==6){
		//only 10 and 01 turned on
		midpoint_1 = vertex_2+(vertex_1-vertex_2)*((DENSITY_TOL-weight_2)/(weight_1-weight_2));
		midpoint_2 = vertex_2+(vertex_4-vertex_2)*((DENSITY_TOL-weight_2)/(weight_4-weight_2));
		tri1 = new Triangle(vertex_2,midpoint_1,midpoint_2);

		midpoint_1 = vertex_3+(vertex_4-vertex_3)*((DENSITY_TOL-weight_3)/(weight_4-weight_3));
		midpoint_2 = vertex_3+(vertex_1-vertex_3)*((DENSITY_TOL-weight_3)/(weight_1-weight_3));
		tri2 = new Triangle(vertex_3,midpoint_1,midpoint_2);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
	}else if(cube_case==7){
		//only 10, 01, and 11 turned on
		midpoint_1 = vertex_2+(vertex_1-vertex_2)*((DENSITY_TOL-weight_2)/(weight_1-weight_2));
		tri1 = new Triangle(vertex_4,vertex_2,midpoint_1);

		midpoint_2 = vertex_3+(vertex_1-vertex_3)*((DENSITY_TOL-weight_3)/(weight_1-weight_3));
		tri2 = new Triangle(vertex_4,midpoint_1,midpoint_2);

		tri3 = new Triangle(vertex_4,midpoint_2,vertex_3);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
		tri_list.push_back(tri3);
	}else if(cube_case==8){
		//only corner 00 is turned on
		tri1 = new Triangle(vertex_1,vertex_3,vertex_2,weight_1,weight_3,weight_2,DENSITY_TOL);
		tri_list.push_back(tri1);
	}else if(cube_case==9){
		//only 00 and 11 turned on
		midpoint_1 = vertex_1+(vertex_3-vertex_1)*((DENSITY_TOL-weight_1)/(weight_3-weight_1));
		midpoint_2 = vertex_1+(vertex_2-vertex_1)*((DENSITY_TOL-weight_1)/(weight_2-weight_1));
		tri1 = new Triangle(vertex_1,midpoint_1,midpoint_2);

		midpoint_1 = vertex_4+(vertex_2-vertex_4)*((DENSITY_TOL-weight_4)/(weight_2-weight_4));
		midpoint_2 = vertex_4+(vertex_3-vertex_4)*((DENSITY_TOL-weight_4)/(weight_3-weight_4));
		tri2 = new Triangle(vertex_4,midpoint_1,midpoint_2);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
	}else if(cube_case==10){
		//only 00 and 10 turned on
		midpoint_1 = vertex_3+(vertex_4-vertex_3)*((DENSITY_TOL-weight_3)/(weight_4-weight_3));
		tri1 = new Triangle(vertex_1,vertex_3,midpoint_1);

		midpoint_2 = vertex_1+(vertex_2-vertex_1)*((DENSITY_TOL-weight_1)/(weight_2-weight_1));
		tri2 = new Triangle(vertex_1,midpoint_1,midpoint_2);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
	}else if(cube_case==11){
		//00, 10, 11
		midpoint_1 = vertex_4+(vertex_2-vertex_4)*((DENSITY_TOL-weight_4)/(weight_2-weight_4));
		tri1 = new Triangle(vertex_3,vertex_4,midpoint_1);

		midpoint_2 = vertex_1+(vertex_2-vertex_1)*((DENSITY_TOL-weight_1)/(weight_2-weight_1));
		tri2 = new Triangle(vertex_3,midpoint_1,midpoint_2);

		tri3 = new Triangle(vertex_3,midpoint_2,vertex_1);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
		tri_list.push_back(tri3);
	}else if(cube_case==12){
		//00 and 01 turned on
		midpoint_1 = vertex_2+(vertex_4-vertex_2)*((DENSITY_TOL-weight_2)/(weight_4-weight_2));
		tri1 = new Triangle(vertex_1,vertex_2,midpoint_1);

		midpoint_2 = vertex_1+(vertex_3-vertex_1)*((DENSITY_TOL-weight_1)/(weight_3-weight_1));
		tri2 = new Triangle(vertex_1,midpoint_1,midpoint_2);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
	}else if(cube_case==13){
		//00, 01, 11
		midpoint_1 = vertex_4+(vertex_3-vertex_4)*((DENSITY_TOL-weight_4)/(weight_3-weight_4));
		tri1 = new Triangle(vertex_2,midpoint_1,vertex_4);

		midpoint_2 = vertex_1+(vertex_3-vertex_1)*((DENSITY_TOL-weight_1)/(weight_3-weight_1));
		tri2 = new Triangle(vertex_2,midpoint_2,midpoint_1);

		tri3 = new Triangle(vertex_2,vertex_1,midpoint_2);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
		tri_list.push_back(tri3);
	}else if(cube_case==14){
		//00, 01, 10
		midpoint_1 = vertex_2+(vertex_4-vertex_2)*((DENSITY_TOL-weight_2)/(weight_4-weight_2));
		tri1 = new Triangle(vertex_1,midpoint_1,vertex_2);

		midpoint_2 = vertex_3+(vertex_4-vertex_3)*((DENSITY_TOL-weight_3)/(weight_4-weight_3));
		tri2 = new Triangle(vertex_1,midpoint_2,midpoint_1);

		tri3 = new Triangle(vertex_1,vertex_3,midpoint_2);

		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
		tri_list.push_back(tri3);
	}else if(cube_case==15){
		tri1 = new Triangle(vertex_1,vertex_3,vertex_2);
		tri2 = new Triangle(vertex_2,vertex_3,vertex_4);
		tri_list.push_back(tri1);
		tri_list.push_back(tri2);
	}
	return tri_list;
}

/*
Run marching cubes algorithm to generate triangles to render. Uses the particles in PARTICLES to
calculate densities at corners of cubes (or squares in 2D). There are 16 cases for squares in
the marching cubes algorithm, 256 for cubes (really 15 up to rotations and reflections).
*/
void marching_cubes(){
	/*we need to break up screen into squares. We should store the squares in an efficient data
	structure so we reuse values previously calculated through iterations.
	*/

	/*look here for cases http://users.polytech.unice.fr/~lingrand/MarchingCubes/algo.html
	For now just making the 16 triangle types as is and enqueuing in list. Would be better to
	make them dynamically.

	Also, only doing uniform marching squares, not doing adaptive squares.
	*/

	VERTEX_MATRIX.clear();
	GRID_DENSITY.clear();
	GRID_BOOL.clear();
	TRIANGLES.clear();

	float error = .0001;
	float step = CUBE_TOL;

	//generate verticies and densities at those verticies
	for (float x = CONTAINER.min.x; x<CONTAINER.max.x+error; x = x+step){
		vector<vector<float> > yz_list;
		vector<vector<Vec3> >vec_yz_list;
		vector<vector<bool> > bool_yz_list;

		for (float y = CONTAINER.min.y; y<CONTAINER.max.y+error; y = y+step){
			vector<float> z_list;
			vector<Vec3> vec_z_list;
			vector<bool> bool_z_list;

			for (float z = CONTAINER.min.z; z<CONTAINER.max.z+error; z = z+step){
				Vec3 vertex(x,y,z);
				float density = density_at_point(vertex);

				z_list.push_back(density);
				vec_z_list.push_back(vertex);
				bool_z_list.push_back(density>DENSITY_TOL);
			}
			yz_list.push_back(z_list);
			vec_yz_list.push_back(vec_z_list);
			bool_yz_list.push_back(bool_z_list);
		}

		GRID_DENSITY.push_back(yz_list);
		GRID_BOOL.push_back(bool_yz_list);
		VERTEX_MATRIX.push_back(vec_yz_list);
	}

	//check marching cubes cases to add new triangles to list. It would be good to use bit operations for speed.
	int m,n,p;//= GRID_DENSITY.size();//this implies uniform. Will need to change this.
	m = n = p = GRID_BOOL.size();

	int i,j,k; //i indicates the row, j indicates the column. 
	i = j = k = 0;

	vector<Triangle*> tri_list;
	Vec3 vertex_1,vertex_2,vertex_3;
	for (int i = 0; i<m-1;i++){
		for (int j = 0; j<n-1;j++){
			for (int k = 0; k<p-1;k++){
				tri_list.clear();
				tri_list = make_triangles(i,j,k);
				TRIANGLES.insert(TRIANGLES.end(),tri_list.begin(),tri_list.end());

				tri_list.clear();
				tri_list = make_triangles(i,j,k+1);
				TRIANGLES.insert(TRIANGLES.end(),tri_list.begin(),tri_list.end());
			}	
		}
	}
	//while (squares_left){ //need to change this later
	//	tri_list.clear();
	//	tri_list = make_triangles(i,j,k);
	//	TRIANGLES.insert(TRIANGLES.end(),tri_list.begin(),tri_list.end());

	//	tri_list.clear();
	//	tri_list = make_triangles(i,j,k+1);
	//	TRIANGLES.insert(TRIANGLES.end(),tri_list.begin(),tri_list.end());

	//	//increment i and j and k or prepare to break loop
	//	if(k==p-2){
	//		if (j==n-2){
	//			if (i==m-2){
	//				squares_left = false;
	//			}else{
	//				j = 0;
	//				i++;
	//			}
	//		}else{
	//			j++;
	//		}
	//	}else{
	//		k++;
	//	}

	//	if (j==n-2){
	//		if (i==m-2){

	//			squares_left = false;
	//		}else{
	//			j = 0;
	//			i++;
	//		}
	//	}else{
	//		j++;
	//	}
	//}
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
