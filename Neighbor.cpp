
#include "Neighbor.h"

#include "Vec3.h"

void Neighbor::reset_box_particles() {
    for (int i = 0; i < box_particles.size(); i++) {
        box_particles[i].clear();
    }
}

int Neighbor::compute_box_num(Vec3 pos) {
    return 9;
}

void Neighbor::place_particles(vector<Particle*> particles, float support_rad) {
    reset_box_particles();
    
    int box_num;
    for (int i = 0; i < particles.size(); i++) {
        // determine box #
        box_num = compute_box_num(particles[i]->position);
        // set box # in particle
        particles[i]->box = box_num;
        // add particle number to corresponding box
        box_particles[box_num].push_back(i);
    }
    
    // get all particles from neighboring boxes and add to particle's neighbor vector
    for (int i = 0; i < particles.size(); i++) {
        particles[i]->neighbors.clear(); // clear out old neighbor vector
        int particle_num;
        
        // need to find out which boxes neighbor the current box
        for (int j = 0; j < box_particles[i].size(); j++) {
            particle_num = box_particles[i][j]; // number of a neighboring particle
            particles[i]->neighbors.push_back(particle_num);
        }
    }
    
}