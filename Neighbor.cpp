
#include "Neighbor.h"

#include <cmath>

#include <iostream>

using namespace std;

void Neighbor::add_to_box_particles(int box_num,int particle_num) {
    box_particles[box_num].push_back(particle_num);
}

void Neighbor::set_particle_neighbors(int particle_num, Particle *p) {
    vector<int> list = box_particles[particle_num];
    for (int i = 0; i < box_particles[particle_num].size(); i++) {
        p->neighbors.push_back(list[i]);
    }
}


int Neighbor::compute_box_num(Vec3 pos, float support_rad, int min, int max) {
    int row = -1,col =-1;
    int width = max - min;
    int blocks_per_row = width / support_rad;
    // find row and column number of box
//    cout<<"Position: "<<pos.x<<", "<<pos.y<<endl;
    
    float curr_x = min, curr_y = min;

	for (int i = 0; i < blocks_per_row && curr_x < float(max); i++) {
		float col_point = abs(curr_x - pos.x);
        float row_point = abs(curr_y - pos.y);
//        cout<<"Col Point: "<<curr_x<<" - "<<pos.x<<" = "<<col_point<<endl;
//        cout<<"Row Point: "<<curr_y<<" - "<<pos.y<<" = "<<row_point<<endl;
		if (col_point <= support_rad) {
			col = i;
		}
        if (row_point <= support_rad) {
			row = i;
		}
        if (row != -1 && col != -1) { break; }
        curr_x += support_rad;
        curr_y += support_rad;
	}
    
    int num = col + row * blocks_per_row;
//    cout<<"Num = "<<num<<" = "<<col<<" + "<<row<<" * "<<blocks_per_row<<endl;
//    cout<<endl;
    if (num > blocks_per_row * blocks_per_row || num < 0) {
        cout<<"Error, incorrect box # assigned in Neighbor: "<<num<<endl;
        num = 0; // temporarily fix bug where y position is ~38
//        exit(0);
    }
    // compute box num given row & col
    return num;
}

void Neighbor::place_particles(vector<Particle*> &particles, float support_rad, Container c) {
    // assuming square container
    int width = c.max.x - c.min.x;
    int min = c.min.x;
    int max = c.max.x;
    box_particles.clear();
    
    for (int i = 0; i < (width/support_rad)*(width/support_rad); i++) {
        box_particles.push_back(vector<int>());
    }
    
    int box_per_row = width/support_rad;
    int box_num;
    for (int i = 0; i < particles.size(); i++) {
        // determine box #
        box_num = compute_box_num(particles[i]->position, support_rad, min, max);
        // set box # in particle
        particles[i]->box = box_num;
        // add particle number to corresponding box
        box_particles[box_num].push_back(i);
    }
    
    // get all particles from neighboring boxes and add to particle's neighbor vector
    for (int i = 0; i < particles.size(); i++) {
        particles[i]->neighbors.clear(); // clear out old neighbor vector
        int particle_num;
        vector<int> neighbor_boxes; // contains the numbers of all the neighboring boxes
        box_num = particles[i]->box;
        
        // each particle is a neighbor to particles in the same box
        neighbor_boxes.push_back(box_num);
        
        if (box_num < box_per_row && box_num != 0 && box_num != box_per_row - 1) {
            // bottom row, not on left or right edge
            neighbor_boxes.push_back(box_num-1);
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row-1);
            neighbor_boxes.push_back(box_num+box_per_row+1);
        } else if (box_num > box_per_row * (box_per_row - 1) - 1 && box_num % box_per_row != 0 && box_num % box_per_row != box_per_row - 1) {
            // top row, not on left or right edge
            neighbor_boxes.push_back(box_num-1);
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row-1);
            neighbor_boxes.push_back(box_num-box_per_row+1);
        } else if (box_num % box_per_row == 0) {
            if (box_num < box_per_row) {
                // bottom left corner
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row-1);
                neighbor_boxes.push_back(box_num+box_per_row+1);
            } else if (box_num > box_per_row * (box_per_row - 1) - 1) {
                // top left corner
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row-1);
                neighbor_boxes.push_back(box_num-box_per_row+1);
            } else {
                // on left side
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row+1);
                neighbor_boxes.push_back(box_num+box_per_row+1);
            }
            
        } else if (box_num % box_per_row == box_per_row - 1) {
            if (box_num < box_per_row) {
                // bottom right corner
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row-1);
            } else if (box_num > box_per_row * (box_per_row - 1) - 1) {
                // top right corner
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row-1);
            } else {
                // on right side
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row-1);
                neighbor_boxes.push_back(box_num+box_per_row-1);
            }
            
        } else {
            // all other locations, not on edges
            neighbor_boxes.push_back(box_num-1);
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row-1);
            neighbor_boxes.push_back(box_num-box_per_row+1);
            neighbor_boxes.push_back(box_num+box_per_row-1);
            neighbor_boxes.push_back(box_num+box_per_row+1);
        }
        
        // need to find out which boxes neighbor the current box
        for (int j = 0; j < neighbor_boxes.size(); j++) {
            int num = neighbor_boxes[j];
            vector<int> neighbor_vec = box_particles[neighbor_boxes[j]]; // particles in a neighboring box
            for (int k = 0; k < neighbor_vec.size(); k++) {
                particle_num = neighbor_vec[k]; // number of a neighboring particle
                Vec3 a = particles[i]->position;
                Vec3 b = particles[k]->position;
                float dist = sqrt(pow((a.x - b.x),2) - pow((a.y - b.y),2));
                if (particle_num != i && dist <= support_rad) {
                    particles[i]->neighbors.push_back(particle_num);
                }
            }
        }
    }
    
}