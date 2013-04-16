
#pragma once
#include "Particle.h"

#include <vector>

#include "Vec3.h"

using namespace std;

#pragma once
class Neighbor {
public:
    vector<vector<int> > box_particles; // vector of particle numbers for each box
    
    Neighbor(void){};
    void place_particles(vector<Particle*>,float);
    void reset_box_particles();
    int compute_box_num(Vec3);
};