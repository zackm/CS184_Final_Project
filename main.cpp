#include "Particle.h"

#include <vector>

using namespace std;

/*******************
* GLOBAL VARIABLES *
*******************/

const float TIMESTEP = .1;
vector<Particle> PARTICLES;
const int NUM_PARTICLES;
const float IDEAL_DENSITY = 1.0f; //for water
const float STIFFNESS = 100.0f; //no idea what it should be set to.

/******************
* Kernel Function *
******************/
float gaussian(Vec3 r_i, Vec3 r_j){
	float distance = r_i.dot(r_j);

	return exp(-4.0f*distance);
}

float gaussian_deriv(Vec3 r_i, Vec3 r_j){
	float distance = r_i.dot(r_j);

	return -8.0f*distance*exp(-4.0f*distance);
}

/*
Create a new list of particles from the old list. Then, throw the old list.
To do this, calculate all quanities in Navier-Stokes, then use timestep to
update particle location from old location and velocity.
*/
void update_particles(){
	vector<Particle> new_particles;

	//update using slow algorithm for now
	Particle base_particle, temp_particle, new_particle;
	Vec3 pressure_gradient;
	float number_density, density, pressure, weight;

	for (int i = 0; i<NUM_PARTICLES; i++){
		density = 0;
		base_particle = PARTICLES[i];

		//for now compute just the term for dv/dt in terms of pressure gradient, ignore the rest.
		for (int j = 0; j<NUM_PARTICLES; j++){
			//update density
			temp_particle = PARTICLES[j];
			number_density += gaussian(base_particle.position,temp_particle.position);
		}
		new_particle.pressure = STIFFNESS*(number_density-IDEAL_DENSITY);
		new_particle.number_density = number_density;


		//for (int j = 0; j<NUM_PARTICLES; j++){

		//	pressure_gradient += density*()*gaussian_derive(
		//}
		//pressure_gradient 
	}

	PARTICLES = new_particles;
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





	return 1;
}
