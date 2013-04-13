#include "Particle.h"

#include <vector>

using namespace std;

/*******************
* GLOBAL VARIABLES *
*******************/

const float TIMESTEP = .1;
vector<Particle> PARTICLES;

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
	vector<Particle> NEW_PARTICLES;

	//update...//

	PARTICLES = NEW_PARTICLES;
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

		alternate form on intel notes for low compressible fluids.
	Step 3:
		pressure gradient at particle j = rho_j*sum(p_i/rho_i^2 + p_j/rho_j^2)*deriv_of_kernel(r_ij)

	Step 4:
		

	*/





	return 1;
}
