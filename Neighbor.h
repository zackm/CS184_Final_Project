
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
    void place_particles(vector<Particle*>,float,Container);
    int compute_box_num(Vec3,float,int,int,int);
};