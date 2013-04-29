
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
    vector<vector<int> > box_particles; // vector of particle numbers for each box
    
    Neighbor(void){};
    void place_particles(vector<Particle*>&,float,Container,int);
    int compute_box_num(Vec3,float,Vec3,Vec3);
    void add_to_box_particles(int,int);
    void set_particle_neighbors(int, Particle*);
};