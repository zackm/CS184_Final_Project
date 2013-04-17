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
vector<vector<float> > GRID_DENSITY;//Grid for marching squares. Probably a better data structure we can use.
vector<vector<bool> > GRID_BOOL; //bools corresponding to that grid
vector<vector<Vec3> > VERTEX_MATRIX;//list of vertices corresponding to the densities on the grid.

const float TIMESTEP = .01;//time elapsed between iterations
const int NUM_PARTICLES = 200;
const Vec3 GRAVITY(0,-9.8f,0);
const float IDEAL_DENSITY = 1000.0f; //for water kg/m^3
const float TEMPERATURE = 293.0f; //kelvin for water at 20 degrees celcius
const float MOLAR_MASS = .0180153f;//for water
const float BOLTZMANN = 8.31446f;//gas constant
const float MASS = 1.0f;//could set it to any number really.
const float STIFFNESS = .01f;//BOLTZMANN*TEMPERATURE/MOLAR_MASS;// for water;
const float VISCOSITY = 1.004f;// for water;
const float COLLISION_RADIUS = .01f;//collision between particles.

const float MAX_KERNEL_RADIUS = .1f;
const float CUBE_TOL = .1f;//either grid size or tolerance for adaptive cubes, reciprocal must be an integer for now.
const float DENSITY_TOL = .5f;//also used for marching grid, for density of the particles

Neighbor NEIGHBOR; //neighbor object used for calculations
const float SUPPORT_RADIUS = .1f;//radius of support used by neighbor function to divide space into grid

bool USE_ADAPTIVE = false; //for adaptive or uniform marching cubes.

const float PI = 3.1415926;

/*
simple dot product between two vectors.
*/
float dot(Vec3 v1, Vec3 v2){
	return v1.x*v2.x+v1.y*v2.y+v1.z+v2.z;
}

/******************	............These derivatives come from wolfram alpha.
* Kernel Function *
******************/
//Poly6 is used for all terms except pressure and viscosity
float poly6_kernel(Vec3 r_i, Vec3 r_j){

	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float h = MAX_KERNEL_RADIUS;

	if (sqrt(mag)>h){
		return 0;
	}

	return (315/(64*PI*pow(h,9.0f)))*pow((h*h - mag),3.0f);
}

Vec3 gradient_kernel(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float h = MAX_KERNEL_RADIUS;

	Vec3 v(0,0,0);
	if (sqrt(mag)>h || sqrt(mag)==0){
		return v;
	}

	float coeff = (45/(PI*pow(h,6.0f)))*pow((h-sqrt(mag)),2.0f)/sqrt(mag);

	return diff_vec*(-coeff);
}

float laplacian_kernel(Vec3 r_i, Vec3 r_j){

	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float h = MAX_KERNEL_RADIUS;

	if (sqrt(mag)>h || sqrt(mag)==0){
		return 0;
	}

	float coeff = 15/(2*PI*pow(h,3.0f));

	return -((6.0f*coeff)/(pow(h,3.0f)*sqrt(mag)))*((-1.0f*h*sqrt(mag))+mag);
}

/********************
* Physics Functions *
********************/
Vec3 kinematic_polynomial(Vec3 acc, Vec3 vel, Vec3 pos,float t){
	return acc*t*t*.5f+vel*t+pos;
}

/*******************
* Particle Methods *
*******************/
void collide_particles(){
	Particle *base_particle, *temp_particle, *new_particle;
	float number_density = 0, density = 0, pressure = 0;

	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];
		vector<int> neighbor_vec = base_particle->neighbors;

	}
}

float density_at_particle(Particle* part){
	float density = poly6_kernel(part->position,part->position);
	Particle* temp_particle;
	vector<int> neighbor_vec = part->neighbors;
	for (int i = 0; i<neighbor_vec.size(); i++){ // changed to neighbors

		temp_particle = PARTICLES[neighbor_vec[i]];
		if (i == neighbor_vec[i]) {
			continue;
		}

		density += temp_particle->mass*poly6_kernel(part->position,temp_particle->position);
	}
	return density;
}

float density_at_point(Vec3 point){
	//first should generate a list of the particles we need to check, using the box for this point.
	int box_number = NEIGHBOR.compute_box_num(point,SUPPORT_RADIUS,CONTAINER.min.x,CONTAINER.max.x);

	Particle* temp_particle;
	vector<int> neighbor_vec = NEIGHBOR.box_particles[box_number];
	float density = 0;
	for (int i = 0; i<neighbor_vec.size(); i++){
		//update density. There will be a problem if the point is exactly equal to sum particle (divide by zero error).
		temp_particle = PARTICLES[neighbor_vec[i]];
		if (i == neighbor_vec[i]) {
			continue;
		}
		density += temp_particle->mass*poly6_kernel(point,temp_particle->position);
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
		density = density_at_particle(base_particle);
		density_list.push_back(density);
		pressure_list.push_back(STIFFNESS*(density-IDEAL_DENSITY));
	}

	//Sets pressure gradient at each point using densities from last loop
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 pressure_gradient(0,0,0);
		
        vector<int> neighbor_vec = base_particle->neighbors;
		int n = neighbor_vec.size();
		for (int j = 0; j<neighbor_vec.size(); j++){ // changed to neighbors
			temp_particle = PARTICLES[neighbor_vec[j]];
            if (i == neighbor_vec[j]) {
                continue;
            }

			Vec3 weight = gradient_kernel(base_particle->position,temp_particle->position);
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

			float weight = laplacian_kernel(base_particle->position,temp_particle->position);
			viscosity_laplacian += ((temp_particle->velocity - base_particle->velocity)/density_list[j])*weight * base_particle->mass;
		}
		viscosity_list.push_back(viscosity_laplacian);
	}

	//Create new particles from old particles and from pressure gradient.
	for (int i = 0; i<NUM_PARTICLES; i++){
		temp_particle = PARTICLES[i];

		//using Navier Stokes, calculate the change in velocity.
		Vec3 acceleration = pressure_grad_list[i]/density_list[i]*(-1.0f)
		+ GRAVITY + (viscosity_list[i]/density_list[i])*VISCOSITY;

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

/************************
* Marching Cube Methods *
************************/
/*
Returns a list of pointers to the different triangles that need to be made for the marching cube case 
corresponding to the input int. Should call some global list of all possible triangles (up to transformation)
to make a particular triangle.
*/
vector<Triangle*> make_triangles(int cube_case, int i, int j){
	//cube define by the corner (i,j), (i,j+1), (i+1,j), and (i+1,j+1).
	Vec3 vertex_1 = VERTEX_MATRIX[i][j];
	Vec3 vertex_2 = VERTEX_MATRIX[i][j+1];
	Vec3 vertex_3 = VERTEX_MATRIX[i+1][j];
	Vec3 vertex_4 = VERTEX_MATRIX[i+1][j+1];
	float weight_1 = GRID_DENSITY[i][j];
	float weight_2 = GRID_DENSITY[i][j+1];
	float weight_3 = GRID_DENSITY[i+1][j];
	float weight_4 = GRID_DENSITY[i+1][j+1];

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

	float error = .001;
	float step = CUBE_TOL;

	//generate verticies and densities at those verticies
	for (float x = CONTAINER.min.x; x<CONTAINER.max.x+error; x = x+step){
		vector<float> y_list;
		vector<Vec3> vec_list;
		vector<bool> bool_list;

		for (float y = CONTAINER.min.y; y<CONTAINER.max.y+error; y = y+step){
			Vec3 vertex(x,y,0);

			float density = density_at_point(vertex);

			y_list.push_back(density);
			vec_list.push_back(vertex);
			bool_list.push_back(density>DENSITY_TOL);
		}

		GRID_DENSITY.push_back(y_list);
		GRID_BOOL.push_back(bool_list);
		VERTEX_MATRIX.push_back(vec_list);
	}

	//check marching cubes cases to add new triangles to list. It would be good to use bit operations for speed.
	const int n = GRID_DENSITY.size();//this implies uniform. Will need to change this.
	const int m = n;

	int i,j; //i indicates the row, j indicates the column. 
	i = j = 0;

	bool squares_left = true;
	vector<Triangle*> tri_list;
	Vec3 vertex_1,vertex_2,vertex_3;
	while (squares_left){ //need to change this later
		tri_list.clear();

		int cube_case = GRID_BOOL[i][j]*8+GRID_BOOL[i][j+1]*4+GRID_BOOL[i+1][j]*2+GRID_BOOL[i+1][j+1];
		tri_list = make_triangles(cube_case,i,j);

		TRIANGLES.insert(TRIANGLES.end(),tri_list.begin(),tri_list.end());

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
		x = float(rand())/(float(RAND_MAX));
		y = float(rand())/(float(RAND_MAX));
		z = 0;//float(rand())/(float(RAND_MAX));

		Vec3 pos(x,y/5.0f,z);
		Vec3 vel(rand()%3,rand()%2,0);
		PARTICLES.push_back(new Particle(pos,vel,MASS));
	}

	glEnable(GL_DEPTH_TEST);

}

void myDisplay(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0,0,400,400);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(90,1.0f,1,-1000);
	glOrtho(CONTAINER.min.x,CONTAINER.max.x,CONTAINER.min.y,CONTAINER.max.y,CONTAINER.min.z,CONTAINER.max.z);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//gluLookAt(0,0,3,0,0,0,0,1,0);
    
	//collide_particles();
	//NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER);

	update_particles();
    NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER);
	//marching_cubes();

	//draw particles
	Particle* temp_part;
	glPointSize(4.0f);
	glBegin(GL_POINTS);
	for (int i = 0; i<PARTICLES.size(); i++){
		temp_part = PARTICLES[i];
		glClearColor(0,0,0,0);
        
        // alternate particle colors depending on box in grid
        if (temp_part->box % 2 == 0) {
            glColor3f(1.0,0,0);
        } else {
            glColor3f(0,1.0,1.0);
        }
		glVertex3f(temp_part->position.x,temp_part->position.y,temp_part->position.z-.1);
	}
	glEnd();

	//draw triangles
	Triangle *temp_triangle;
	glPointSize(4.0f);
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














//Vec3 poly6_gradient(Vec3 r_i, Vec3 r_j){
//
//	Vec3 diff_vec = r_i-r_j;
//	float mag = dot(diff_vec,diff_vec);
//
//	float h = MAX_KERNEL_RADIUS;
//
//	Vec3 v(0,0,0);
//	if (sqrt(mag)>h){
//		return v;
//	}
//
//
//}
//
//float poly6_laplacian(Vec3 r_i, Vec3 r_j){
//}

//spiky is used for pressure
//float spiky_kernel(Vec3 r_i, Vec3 r_j){
//
//	Vec3 diff_vec = r_i-r_j;
//	float mag = dot(diff_vec,diff_vec);
//
//	float h = MAX_KERNEL_RADIUS;
//
//	if (sqrt(mag)>h){
//		return 0;
//	}
//
//	return (15/(PI*pow(h,6.0f)))*pow((h-sqrt(mag)),3.0f);
//}

//float spiky_laplacian(Vec3 r_i, Vec3 r_j){
//	Vec3 diff_vec = r_i-r_j;
//	float mag = dot(diff_vec,diff_vec);
//
//	float h = MAX_KERNEL_RADIUS;
//
//	if (sqrt(mag)>h){
//		return 0;
//	}
//
//	return (90/(PI*pow(h,6.0f)))*((-h*h/sqrt(mag))+3.0f*h - 2.0f*sqrt(mag));
//}

////viscosity kernel
//float viscosity_kernel(Vec3 r_i, Vec3 r_j){
//
//	Vec3 diff_vec = r_i-r_j;
//	float mag = dot(diff_vec,diff_vec);
//
//	float h = MAX_KERNEL_RADIUS;
//
//	if (sqrt(mag)>h){
//		return 0;
//	}
//
//	return (15/(2*PI*h*h*h))*(-.5f*(pow(mag/h,3.0f))+(pow(mag/h,2.0f))+.5f*(h/mag)-1);
//}
//
//Vec3 viscosity_gradient(Vec3 r_i, Vec3 r_j){
//}

//Gaussian is an extra kernel type we could use.
//float gaussian(Vec3 r_i, Vec3 r_j){
//	//should add another parameter (max distance value)
//	Vec3 diff_vec = r_i-r_j;
//	float mag = dot(diff_vec,diff_vec);
//
//	if (sqrt(mag)>MAX_KERNEL_RADIUS){
//		return 0;
//	}
//
//	return exp(-4.0f*mag)/sqrt(2*PI);
//}
//
//Vec3 gaussian_grad(Vec3 r_i, Vec3 r_j){
//	Vec3 diff_vec = r_i-r_j;
//	float mag = dot(diff_vec,diff_vec);
//
//	Vec3 v(0,0,0);
//	if (sqrt(mag)>MAX_KERNEL_RADIUS){
//		return v;
//	}
//
//	float coeff = -8.0f*exp(-4.0f*mag);
//
//	Vec3 grad(diff_vec.x,diff_vec.y,diff_vec.z);
//
//	return coeff*grad/sqrt(2*PI);
//}
//
//float gaussian_laplacian(Vec3 r_i, Vec3 r_j){
//	Vec3 diff_vec = r_i-r_j;
//	float mag = dot(diff_vec,diff_vec);
//
//	if (sqrt(mag)>MAX_KERNEL_RADIUS){
//		return 0;
//	}
//
//	float coeff = 8.0f*exp(-4.0f*mag);
//
//	//I got this formula from wolfram alpha, so may want to double check it later.
//	return coeff * (8.0*(mag)-3)/sqrt(2*PI);
//}


		//vector<float> first_row_density = GRID_DENSITY[i];
		//vector<Vec3> first_row_vertices = VERTEX_MATRIX[i];
		//vector<float> second_row_density = GRID_DENSITY[i+1];
		//vector<Vec3> second_row_vertices = VERTEX_MATRIX[i+1];

		//density[0][0] = 
		//float density_00 = first_row_density[j];
		//float density_01 = first_row_density[j+1];
		//float density_10 = second_row_density[j];
		//float GRID_DENSITY[i+1][j+1] = second_row_density[j+1];

		//bool corner_00 = (GRID_DENSITY[i][j]>DENSITY_TOL),
		//	corner_01 = (GRID_DENSITY[i][j+1]>DENSITY_TOL),
		//	GRID_BOOL[i+1][j] = (GRID_DENSITY[i+1][j]>DENSITY_TOL),
		//	corner_11 = (GRID_DENSITY[i+1][j+1]>DENSITY_TOL);

		/*int corner_0 = (density_00>DENSITY_TOL) ? 1 : 0,
			corner_1 = (density_01>DENSITY_TOL) ? 1 : 0,
			corner_2 = (density_10>DENSITY_TOL) ? 1 : 0,
			corner_3 = (GRID_DENSITY[i+1][j+1]>DENSITY_TOL) ? 1 : 0;
			*/
		////now check the case and add the corresponding triangle type.
		//if (GRID_BOOL[0][0]){
		//	if (GRID_BOOL[i][j+1]){
		//		if (GRID_BOOL[i+1][j]){
		//			if (GRID_BOOL[i+1][j+1]){
		//				//make whole square out of two triangles.
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i][j],VERTEX_MATRIX[i+1][j],VERTEX_MATRIX[i][j+1]);
		//				Triangle* tri2 = new Triangle(VERTEX_MATRIX[i+1][j],VERTEX_MATRIX[i+1][j+1],VERTEX_MATRIX[i][j+1]);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//			}else{
		//				vertex_1 = VERTEX_MATRIX[i][j];
		//				vertex_2 = VERTEX_MATRIX[i+1][j];
		//				vertex_3 = VERTEX_MATRIX[i+1][j+1];
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i][j],VERTEX_MATRIX[i+1][j],(VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i+1][j+1])/2.0f);
		//				
		//				vertex_2 = ;
		//				vertex_3 = ;
		//				Triangle* tri2 = new Triangle((VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i+1][j+1])/2.0f,(VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i+1][j+1])/2.0f,VERTEX_MATRIX[i][j]);
		//				
		//				vertex_2 = ;
		//				vertex_3 = ;
		//				Triangle* tri3 = new Triangle((VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i+1][j+1])/2.0f,VERTEX_MATRIX[i][j+1],VERTEX_MATRIX[i][j]);
		//				
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//				TRIANGLES.push_back(tri3);
		//			}
		//		}else{
		//			if (GRID_BOOL[i+1][j+1]){
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i][j+1],VERTEX_MATRIX[i][j],(VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i][j])/2.0f);
		//				Triangle* tri2 = new Triangle((VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i][j])/2.0f,(VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i+1][j+1])/2.0f,VERTEX_MATRIX[i][j+1]);
		//				Triangle* tri3 = new Triangle((VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i+1][j+1])/2.0f,VERTEX_MATRIX[i+1][j+1],VERTEX_MATRIX[i][j+1]);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//				TRIANGLES.push_back(tri3);
		//			}else{
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i][j],(VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i+1][j])/2,VERTEX_MATRIX[i][j+1]);
		//				Triangle* tri2 = new Triangle(VERTEX_MATRIX[i][j+1],(VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i+1][j])/2,(VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i+1][j+1])/2.0f);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//			}
		//		}
		//	}else{
		//		if (GRID_BOOL[i+1][j]){
		//			if (GRID_BOOL[i+1][j+1]){
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i+1][j],VERTEX_MATRIX[i+1][j+1],(VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i+1][j+1])/2.0f);
		//				Triangle* tri2 = new Triangle((VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i+1][j+1])/2.0f,(VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i][j+1])/2.0f,VERTEX_MATRIX[i+1][j]);
		//				Triangle* tri3 = new Triangle((VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i][j+1])/2.0f,VERTEX_MATRIX[i][j],VERTEX_MATRIX[i+1][j]);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//				TRIANGLES.push_back(tri3);
		//			}else{
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i][j],VERTEX_MATRIX[i+1][j],(VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i+1][j+1])/2.0f);
		//				Triangle* tri2 = new Triangle((VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i+1][j+1])/2.0f,(VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i][j+1])/2.0f,VERTEX_MATRIX[i][j]);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//			}
		//		}else{
		//			if (GRID_BOOL[i+1][j+1]){
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i][j],(VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i+1][j])/2.0f,(VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i][j+1])/2.0f);
		//				Triangle* tri2 = new Triangle(VERTEX_MATRIX[i+1][j+1],(VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i+1][j+1])/2.0f,(VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i+1][j+1])/2.0f);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//			}else{
		//				vertex_1 = VERTEX_MATRIX[i][j];
		//				vertex_2 = VERTEX_MATRIX[i+1][j];
		//				vertex_3 = VERTEX_MATRIX[i][j+1];
		//				Triangle* tri = new Triangle(vertex_1,vertex_2,vertex_3,GRID_DENSITY[i][j],GRID_DENSITY[i+1][j],GRID_DENSITY[i][j+1],DENSITY_TOL);
		//				TRIANGLES.push_back(tri);
		//			}
		//		}
		//	}
		//}else{
		//	if (GRID_BOOL[i][j+1]){
		//		if (GRID_BOOL[i+1][j]){
		//			if (GRID_BOOL[i+1][j+1]){
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i+1][j+1],VERTEX_MATRIX[i][j+1],(VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i][j])/2.0f);
		//				Triangle* tri2 = new Triangle((VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i][j])/2.0f,(VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i+1][j])/2.0f,VERTEX_MATRIX[i+1][j+1]);
		//				Triangle* tri3 = new Triangle((VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i+1][j])/2.0f,VERTEX_MATRIX[i+1][j],VERTEX_MATRIX[i+1][j+1]);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//				TRIANGLES.push_back(tri3);
		//			}else{
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i][j+1],(VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i][j+1])/2.0f,(VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i+1][j+1])/2.0f);
		//				Triangle* tri2 = new Triangle(VERTEX_MATRIX[i+1][j],(VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i+1][j+1])/2.0f,(VERTEX_MATRIX[i][j]+VERTEX_MATRIX[i+1][j])/2.0f);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//			}
		//		}else{
		//			if (GRID_BOOL[i+1][j+1]){
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i][j+1],(VERTEX_MATRIX[i][j+1]+VERTEX_MATRIX[i][j])/2.0f,(VERTEX_MATRIX[i+1][j+1]+VERTEX_MATRIX[i+1][j])/2.0f);
		//				Triangle* tri2 = new Triangle((VERTEX_MATRIX[i+1][j+1]+VERTEX_MATRIX[i+1][j])/2.0f, VERTEX_MATRIX[i+1][j+1],VERTEX_MATRIX[i][j+1]);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//			}else{
		//				vertex_1 = VERTEX_MATRIX[i][j+1];
		//				vertex_2 = VERTEX_MATRIX[i][j];
		//				vertex_3 = VERTEX_MATRIX[i+1][j+1];
		//				Triangle* tri = new Triangle(vertex_1,vertex_2,vertex_3,GRID_DENSITY[i][j+1],GRID_DENSITY[i][j],GRID_DENSITY[i+1][j+1],DENSITY_TOL);
		//				
		//				TRIANGLES.push_back(tri);
		//			}
		//		}
		//	}else{
		//		if (GRID_BOOL[i+1][j]){
		//			if (GRID_BOOL[i+1][j+1]){
		//				Triangle* tri1 = new Triangle(VERTEX_MATRIX[i+1][j],VERTEX_MATRIX[i+1][j+1],(VERTEX_MATRIX[i+1][j+1]+VERTEX_MATRIX[i][j+1])/2.0f);
		//				Triangle* tri2 = new Triangle((VERTEX_MATRIX[i+1][j+1]+VERTEX_MATRIX[i][j+1])/2.0f,(VERTEX_MATRIX[i+1][j]+VERTEX_MATRIX[i][j])/2.0f,VERTEX_MATRIX[i+1][j]);
		//				TRIANGLES.push_back(tri1);
		//				TRIANGLES.push_back(tri2);
		//			}else{
		//				vertex_1 = VERTEX_MATRIX[i+1][j];
		//				vertex_2 = VERTEX_MATRIX[i+1][j+1];
		//				vertex_3 = VERTEX_MATRIX[i][j];
		//				Triangle* tri = new Triangle(vertex_1,vertex_2,vertex_3,GRID_DENSITY[i+1][j],GRID_DENSITY[i+1][j+1],GRID_DENSITY[i][j],DENSITY_TOL);
		//				TRIANGLES.push_back(tri);
		//			}
		//		}else{
		//			if (GRID_BOOL[i+1][j+1]){
		//				vertex_1 = VERTEX_MATRIX[i+1][j+1];
		//				vertex_2 = VERTEX_MATRIX[i][j+1];
		//				vertex_3 = VERTEX_MATRIX[i+1][j];
		//				Triangle* tri = new Triangle(vertex_1,vertex_2,vertex_3,GRID_DENSITY[i+1][j+1],GRID_DENSITY[i][j+1],GRID_DENSITY[i+1][j],DENSITY_TOL);
		//				TRIANGLES.push_back(tri);
		//			}else{
		//				//do nothing, no corners are turned on.
		//			}
		//		}
		//	}
		//}


///************
//* Overloads *
//************/
//Vec3 operator * (float t, const Vec3& arg_vec){
//	Vec3 new_vec;
//	new_vec.x = arg_vec.x*t;
//	new_vec.y = arg_vec.y*t;
//	new_vec.z = arg_vec.z*t;
//	return new_vec;
//}
//
////Vec3 operator * (const Vec3& arg_vec,float t){
////	Vec3 new_vec;
////	new_vec.x = arg_vec.x*t;
////	new_vec.y = arg_vec.y*t;
////	new_vec.z = arg_vec.z*t;
////	return new_vec;
////}
//
//Vec3 operator / (float t, const Vec3& arg_vec){
//	Vec3 new_vec;
//	new_vec.x = arg_vec.x/t;
//	new_vec.y = arg_vec.y/t;
//	new_vec.z = arg_vec.z/t;
//	return new_vec;
//}
//
//Vec3 operator / (const Vec3& arg_vec,float t){
//	Vec3 new_vec;
//	new_vec.x = arg_vec.x/t;
//	new_vec.y = arg_vec.y/t;
//	new_vec.z = arg_vec.z/t;
//	return new_vec;
//}
