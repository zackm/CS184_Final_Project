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

Container CONTAINER(Vec3(2,2,2),Vec3(-2,-2,-2));//very simple cube for now. Later, make it a particle itself.
vector<Particle*> PARTICLES;//particles that we do SPH on.
vector<Triangle*> TRIANGLES;//triangles made from marching cubes to render
vector<vector<float> > GRID_DENSITY;//Grid for marching squares. Probably a better data structure we can use.
vector<vector<Vec3> > VERTEX_MATRIX;//list of vertices corresponding to the densities on the grid.
const float TIMESTEP = .05;//time elapsed between iterations
const int NUM_PARTICLES = 20;
const Vec3 GRAVITY(0,-9.8f,0);
const float IDEAL_DENSITY = 1.0f; //for water
const float STIFFNESS = 1.0f; //no idea what it should be set to for water.
const float VISCOSITY = .01f;

const float CUBE_TOL = .1f;//either grid size or tolerance for adaptive cubes.
const float DENSITY_TOL = 1.0f;//also used for marching grid, for density of the particles

Neighbor NEIGHBOR; //neighbor object used for calculations
const float SUPPORT_RADIUS = 10.0f;//radius of support used by neighbor function to divide space into grid

bool USE_ADAPTIVE = false; //for adaptive or uniform marching cubes.

/************
* Overloads *
************/
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
//need to add coefficient so that we normalize kernel.
float gaussian(Vec3 r_i, Vec3 r_j){
	//should add another parameter (max distance value)
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	return exp(-4.0f*mag);
}

Vec3 gaussian_grad(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float coeff = -8.0f*exp(-4.0f*mag);

	Vec3 grad(diff_vec.x,diff_vec.y,diff_vec.z);

	return coeff*grad;
}

float gaussian_laplacian(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float coeff = 8.0f*exp(-4.0f*mag);

	//I got this formula from wolfram alpha, so may want to double check it later.
	return coeff * (8.0*(mag)-3);
}

/********************
* Physics Functions *
********************/
Vec3 kinematic_polynomial(Vec3 acc, Vec3 vel, Vec3 pos,float t){
	return .5f*acc*t*t+vel*t+pos;
}

float density_at_point(Vec3 point){
	float density = 0;
	for (int i = 0; i<NUM_PARTICLES; i++){
		//update density. There will be a problem if the point is exactly equal to sum particle (divide by zero error).
		Particle *temp_particle = PARTICLES[i];

		density += temp_particle->mass*gaussian(point,temp_particle->position);
	}

	return density;
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
	vector<Vec3> viscosity_list;

	//update using slow algorithm for now
	Particle *base_particle, *temp_particle, *new_particle;
	float number_density = 0, density = 0, pressure = 0;

	//Sets density at each point
	for (int i = 0; i<NUM_PARTICLES; i++){
		density = 0;
		base_particle = PARTICLES[i];
		density += density_at_point(base_particle->position);
		density_list.push_back(density);
		pressure_list.push_back(STIFFNESS*(density-IDEAL_DENSITY));
	}

	//Sets pressure gradient at each point using densities from last loop
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 pressure_gradient(0,0,0);
		for (int j = 0; j<NUM_PARTICLES && j!=i; j++){ // change to neighbors
			temp_particle = PARTICLES[j];

			Vec3 weight = gaussian_grad(base_particle->position,temp_particle->position);
			pressure_gradient += temp_particle->mass * ((pressure_list[i]+pressure_list[j])/(2.0f*density_list[j]))*weight; 

		}
		pressure_grad_list.push_back(pressure_gradient);
	}

	//Sets viscosity at each particle
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 viscosity_laplacian(0,0,0);
		for (int j = 0; j<NUM_PARTICLES; j++){ // change to neighbors
			temp_particle = PARTICLES[j];

			float weight = gaussian_laplacian(base_particle->position,temp_particle->position);
			viscosity_laplacian += base_particle->mass*((base_particle->velocity - temp_particle->velocity)/density_list[j])*weight;
		}
		viscosity_list.push_back(viscosity_laplacian);
	}

	//Create new particles from old particles and from pressure gradient.
	for (int i = 0; i<NUM_PARTICLES; i++){
		temp_particle = PARTICLES[i];

		//using Navier Stokes, calculate the change in velocity.
		Vec3 acceleration = -1.0f*pressure_grad_list[i]/density_list[i]
		+ GRAVITY + VISCOSITY * viscosity_list[i]/density_list[i];

		Vec3 velocity = temp_particle->velocity;
		Vec3 position = temp_particle->position;
		Vec3 new_position = kinematic_polynomial(acceleration,velocity,position,TIMESTEP);
		Vec3 new_velocity = temp_particle->velocity + acceleration*TIMESTEP;
		float mass = temp_particle->mass;
		temp_particle = new Particle(new_position,new_velocity,mass);

		CONTAINER.in_container(temp_particle); //applies reflections if outside of boundary.

		new_particles.push_back(temp_particle);
	}
	PARTICLES = new_particles;
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
	TRIANGLES.clear();

	float error = .001;
	float step = CUBE_TOL;
	//float min_x = CONTAINER.min.x;
	//float min_y = CONTAINER.min.y;
	//float max_x = CONTAINER.max.x;
	//float max_y = CONTAINER.max.y;

	//generate verticies and densities at those verticies
	for (float x = CONTAINER.min.x; x<CONTAINER.max.x+error; x = x+step){
		vector<float> y_list;
		vector<Vec3> vec_list;

		for (float y = CONTAINER.min.y; y<CONTAINER.max.y+error; y = y+step){
			Vec3 vertex(x,y,0);

			float density = density_at_point(vertex);

			y_list.push_back(density);
			vec_list.push_back(vertex);
		}

		GRID_DENSITY.push_back(y_list);
		VERTEX_MATRIX.push_back(vec_list);
	}

	//check marching cubes cases to add new triangles to list. It would be good to use bit operations for speed.
	int n = GRID_DENSITY.size();//this implies uniform. Will need to change this.
	int m = n;

	int i,j; //i indicates the row, j indicates the column. 
	i = j = 0;

	bool squares_left = true;
	while (squares_left){ //need to change this later

		vector<float> first_row_density = GRID_DENSITY[i];
		vector<Vec3> first_row_vertices = VERTEX_MATRIX[i];
		vector<float> second_row_density = GRID_DENSITY[i+1];
		vector<Vec3> second_row_vertices = VERTEX_MATRIX[i+1];

		float density_00 = first_row_density[j];
		float density_01 = first_row_density[j+1];
		float density_10 = second_row_density[j];
		float density_11 = second_row_density[j+1];

		bool corner_00 = density_00>DENSITY_TOL,
			 corner_01 = density_01>DENSITY_TOL,
			 corner_10 = density_10>DENSITY_TOL,
			 corner_11 = density_11>DENSITY_TOL;

		//now check the case and add the corresponding triangle type.
		if (corner_00){
			if (corner_01){
				if (corner_10){
					if (corner_11){
						//make whole square out of two triangles.
						Triangle* tri1 = new Triangle(first_row_vertices[j],second_row_vertices[j],first_row_vertices[j+1]);
						Triangle* tri2 = new Triangle(second_row_vertices[j],second_row_vertices[j+1],first_row_vertices[j+1]);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
					}else{
						Triangle* tri1 = new Triangle(first_row_vertices[j],second_row_vertices[j],(second_row_vertices[j]+second_row_vertices[j+1])/2.0f);
						Triangle* tri2 = new Triangle((second_row_vertices[j]+second_row_vertices[j+1])/2.0f,(first_row_vertices[j+1]+second_row_vertices[j+1])/2.0f,first_row_vertices[j]);
						Triangle* tri3 = new Triangle((first_row_vertices[j+1]+second_row_vertices[j+1])/2.0f,first_row_vertices[j+1],first_row_vertices[j]);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
						TRIANGLES.push_back(tri3);
					}
				}else{
					if (corner_11){
						Triangle* tri1 = new Triangle(first_row_vertices[j+1],first_row_vertices[j],(second_row_vertices[j]+first_row_vertices[j])/2.0f);
						Triangle* tri2 = new Triangle((second_row_vertices[j]+first_row_vertices[j])/2.0f,(second_row_vertices[j]+second_row_vertices[j+1])/2.0f,first_row_vertices[j+1]);
						Triangle* tri3 = new Triangle((second_row_vertices[j]+second_row_vertices[j+1])/2.0f,second_row_vertices[j+1],first_row_vertices[j+1]);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
						TRIANGLES.push_back(tri3);
					}else{
						Triangle* tri1 = new Triangle(first_row_vertices[j],(first_row_vertices[j]+second_row_vertices[j])/2,first_row_vertices[j+1]);
						Triangle* tri2 = new Triangle(second_row_vertices[j+1],(first_row_vertices[j]+second_row_vertices[j])/2,second_row_vertices[j+1]);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
					}
				}
			}else{
				if (corner_10){
					if (corner_11){
						Triangle* tri1 = new Triangle(second_row_vertices[j],second_row_vertices[j+1],(first_row_vertices[j+1]+second_row_vertices[j+1])/2.0f);
						Triangle* tri2 = new Triangle((first_row_vertices[j+1]+second_row_vertices[j+1])/2.0f,(first_row_vertices[j]+first_row_vertices[j+1])/2.0f,second_row_vertices[j]);
						Triangle* tri3 = new Triangle((first_row_vertices[j]+first_row_vertices[j+1])/2.0f,first_row_vertices[j],second_row_vertices[j]);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
						TRIANGLES.push_back(tri3);
					}else{
						Triangle* tri1 = new Triangle(first_row_vertices[j],second_row_vertices[j],(second_row_vertices[j]+second_row_vertices[j+1])/2.0f);
						Triangle* tri2 = new Triangle((second_row_vertices[j]+second_row_vertices[j+1])/2.0f,(first_row_vertices[j]+first_row_vertices[j+1])/2.0f,first_row_vertices[j]);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
					}
				}else{
					if (corner_11){
						Triangle* tri1 = new Triangle(first_row_vertices[j],(first_row_vertices[j]+second_row_vertices[j])/2.0f,(first_row_vertices[j]+first_row_vertices[j+1])/2.0f);
						Triangle* tri2 = new Triangle(second_row_vertices[j+1],(first_row_vertices[j+1]+second_row_vertices[j+1])/2.0f,(second_row_vertices[j]+second_row_vertices[j+1])/2.0f);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
					}else{
						Triangle* tri = new Triangle(first_row_vertices[j],(second_row_vertices[j]+first_row_vertices[j])/2.0f,(first_row_vertices[j+1]+first_row_vertices[j])/2.0f);
						TRIANGLES.push_back(tri);
					}
				}
			}
		}else{
			if (corner_01){
				if (corner_10){
					if (corner_11){
						Triangle* tri1 = new Triangle(second_row_vertices[j+1],first_row_vertices[j+1],(first_row_vertices[j+1]+first_row_vertices[j])/2.0f);
						Triangle* tri2 = new Triangle((first_row_vertices[j+1]+first_row_vertices[j])/2.0f,(first_row_vertices[j]+second_row_vertices[j])/2.0f,second_row_vertices[j+1]);
						Triangle* tri3 = new Triangle((first_row_vertices[j]+second_row_vertices[j])/2.0f,second_row_vertices[j],second_row_vertices[j+1]);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
						TRIANGLES.push_back(tri3);
					}else{
						Triangle* tri1 = new Triangle(first_row_vertices[j+1],(first_row_vertices[j]+first_row_vertices[j+1])/2.0f,(first_row_vertices[j+1]+second_row_vertices[j+1])/2.0f);
						Triangle* tri2 = new Triangle(second_row_vertices[j],(second_row_vertices[j]+second_row_vertices[j+1])/2.0f,(first_row_vertices[j]+second_row_vertices[j])/2.0f);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
					}
				}else{
					if (corner_11){
						Triangle* tri1 = new Triangle(first_row_vertices[j+1],(first_row_vertices[j+1]+first_row_vertices[j])/2.0f,(second_row_vertices[j+1]+second_row_vertices[j])/2.0f);
						Triangle* tri2 = new Triangle((second_row_vertices[j+1]+second_row_vertices[j])/2.0f, second_row_vertices[j+1],first_row_vertices[j+1]);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
					}else{
						Triangle* tri = new Triangle((first_row_vertices[j]+first_row_vertices[j+1])/2.0f,(second_row_vertices[j+1]+first_row_vertices[j+1])/2.0f,first_row_vertices[j+1]);
						TRIANGLES.push_back(tri);
					}
				}
			}else{
				if (corner_10){
					if (corner_11){
						Triangle* tri1 = new Triangle(second_row_vertices[j],second_row_vertices[j+1],(second_row_vertices[j+1]+first_row_vertices[j+1])/2.0f);
						Triangle* tri2 = new Triangle((second_row_vertices[j+1]+first_row_vertices[j+1])/2.0f,(second_row_vertices[j]+first_row_vertices[j])/2.0f,second_row_vertices[j]);
						TRIANGLES.push_back(tri1);
						TRIANGLES.push_back(tri2);
					}else{
						Triangle* tri = new Triangle((first_row_vertices[j]+second_row_vertices[j])/2.0f,second_row_vertices[j],(second_row_vertices[j+1]+second_row_vertices[j])/2.0f);
						TRIANGLES.push_back(tri);
					}
				}else{
					if (corner_11){
						Triangle* tri = new Triangle(second_row_vertices[j],(second_row_vertices[j+1]+second_row_vertices[j])/2.0f,(first_row_vertices[j+1]+second_row_vertices[j])/2.0f);
						TRIANGLES.push_back(tri);
					}else{
						//do nothing, no corners are turned on.
					}
				}
			}
		}


		//increment i and j or prepare to break loop
		if (j==n-2){
			if (i==m-2){
				squares_left = false;
			}else{
				j = 0;
				i++;
			}
		}else{
			j++;
		}
	}
}

/*****************
* OpenGL Methods *
*****************/
void initScene(){
	//create a list of random particles
	float x,y,z;
	for (int i = 0; i<NUM_PARTICLES; i++){
		x = float(rand())/(float(RAND_MAX))-1;
		y = float(rand())/(float(RAND_MAX))-1;
		z = 0;//float(rand())/(float(RAND_MAX));
		Vec3 pos(x,y,z);
		Vec3 vel(rand()%3,rand()%3,0);
		float mass = 1.0f;
		PARTICLES.push_back(new Particle(pos,vel,mass));
	}


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
}

void myDisplay(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0,0,400,400);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(90,1.0f,1,-1000);
	glOrtho(-2,2,-2,2,1,-1);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//gluLookAt(0,0,3,0,0,0,0,1,0);

    //NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS);
	update_particles();
	marching_cubes();

	Triangle *temp_triangle;
	glPointSize(4.0f);
	for (int i = 0; i<TRIANGLES.size(); i++){
		temp_triangle = TRIANGLES[i];
		glClearColor(0,0,0,0);
		glColor3f(0,0,1.0);
		glBegin(GL_TRIANGLES);
		glVertex3f(temp_triangle->a.x,temp_triangle->a.y,temp_triangle->a.z);
		glVertex3f(temp_triangle->b.x,temp_triangle->b.y,temp_triangle->b.z);
		glVertex3f(temp_triangle->c.x,temp_triangle->c.y,temp_triangle->c.z);
		glEnd();
	}
	glPopMatrix();

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
	glutCreateWindow("Tyler and Zack Final Project");

	initScene();

	glutDisplayFunc(myDisplay);
	glutIdleFunc(myDisplay);


	glutMainLoop();

	return 1;
}
