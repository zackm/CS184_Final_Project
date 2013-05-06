
#pragma once
#include "Particle.h"

#pragma once
#include "Container.h"

#include <vector>

#pragma once
#include "Vec3.h"

using namespace std;

#pragma once
class Neighbor {
public:
    vector<vector<int> > box_particles; // stores vector of particle numbers that are containted in a specific box number
    vector<vector<int> > box_neighbors; // stores vector of particle numbers that neighbor a specific box number (including particles inside the box)
    
    Neighbor(void){};
    void place_particles(vector<Particle*>&,float,Container,int,bool);
    int compute_box_num(Vec3,float,float,float);
	int compute_box_num(Vec3,float,float,float,bool);
    void add_to_box_neighbors(int,int);
    void set_particle_neighbors(int, Particle*);
};