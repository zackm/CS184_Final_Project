
#include "Neighbor.h"

#include <cmath>

#include <iostream>

#include <algorithm>

using namespace std;

void Neighbor::add_to_box_neighbors(int box_num,int particle_num) {
    bool contains = false;
    vector<int> current_neighbors = box_neighbors[box_num];
    for (int i = 0; i < current_neighbors.size(); i++) {
        if (current_neighbors[i] == particle_num) {
            contains = true;
            break; // the current particle is already in the neighbor vector for this box number
        }
    }
	if (!contains) {
        box_neighbors[box_num].push_back(particle_num);
    }
}

void Neighbor::set_particle_neighbors(int particle_num, Particle *p) {
	vector<int> list = box_particles[particle_num];
    
	for (int i = 0; i < box_particles[particle_num].size(); i++) {
		p->neighbors.push_back(list[i]);
	}
}

float dot(Vec3,Vec3);

/*
 Compute Box Number Function:
 Given a particle position, support radius, maximum container x value, and minimum container x value, returns
 the box number for the neighbor algorithm. Assumes the container is a cube and inverse of support radius is
 an integer.
 */
int Neighbor::compute_box_num(Vec3 pos, float support_rad, float max_point, float min_point) {
	int col = floor(pos.x/support_rad);
	int row = floor(pos.y/support_rad);
	int depth = floor(pos.z/support_rad);

	//int row = -1,col = -1, depth = -1;
    float width = max_point - min_point;
    int box_per_row = (int)(width / support_rad); // casting to int, assuming support radius evenly divides width
    
 //   // current x,y,z locations of cell traversal
 //   float curr_x = min_point, curr_y = min_point, curr_z = min_point;
 //   
 //   // holds the difference between current x,y,z traversal locations and the particle x,y,z locations
	//float col_point, row_point, depth_point;
 //   
 //   // loop until assign box number (1D numbering) in x,y, and z directions
	//for (int i = 0; i < box_per_row && curr_x < max_point; i++) {
 //       col_point = abs(pos.x - curr_x-support_rad/2.0f);
 //       row_point = abs(pos.y - curr_y-support_rad/2.0f);
 //       depth_point = abs(pos.z - curr_z-support_rad/2.0f);
 //       
	//	if (col_point <= support_rad && col == -1) {
	//		col = i;
	//	}
 //       if (row_point <= support_rad && row == -1) {
	//		row = i;
	//	}
 //       if (depth_point <= support_rad && depth == -1) {
 //           depth = i;
 //       }
 //       
 //       if (row != -1 && col != -1 && depth != -1) {
	//		break;
	//	}
 //       curr_x += support_rad;
 //       curr_y += support_rad;
 //       curr_z += support_rad;
	//}
 //   

    // combine box numbers into 3D numbering
    int num = row + col * box_per_row + depth * box_per_row * box_per_row;
    
    if (num >= box_per_row * box_per_row * box_per_row || num < 0 || row == -1 || col == -1 || depth == -1) {
        //cout<<"Error, incorrect box # assigned in Neighbor: "<<num<<endl;
        num = 0; // set box number to 0 to prevent bad vector access
    }
    return num;
}

/*
 Compute Box Number Function (for raytracer)
 This is the same as the other compute box number function, except that it handles errors differently. If
 a particle position is outside the min and max point range, it return -1. Same assumptions hold as before.
 */
int Neighbor::compute_box_num(Vec3 pos, float support_rad, float max_point, float min_point, bool flag) {
//    int row = -1,col = -1, depth = -1;
    int col = floor(pos.x/support_rad);
	int row = floor(pos.y/support_rad);
	int depth = floor(pos.z/support_rad);

	//int row = -1,col = -1, depth = -1;
    float width = max_point - min_point;
    int box_per_row = (int)(width / support_rad);
    
    int num = col + row * box_per_row + depth * box_per_row * box_per_row;
    
    if (num >= box_per_row * box_per_row * box_per_row || num < 0 || row == -1 || col == -1 || depth == -1) {
        //cout<<"Error, incorrect box # assigned in Neighbor: "<<num<<endl;
        num = -1; // set box number to -1 to signify out of container
        // exit(0);
    }
    //compute box num given row & col
    return num;
}
/*
 Place Particles Function:
 Given a vector of particles, a support radius, a container, and the number of particles, this will assign a box
 number to each particle, along with a vector of neighboring particle numbers. Works by first 
 */
void Neighbor::place_particles(vector<Particle*>& particles,float support_rad, Container c, int num_particles, bool SURFACE){
    float width = c.max.x - c.min.x;
    float min = c.min.x;
    float max = c.max.x;
    int box_per_row = (int)(width/support_rad); // numbers of boxes along each axis
    int square_face = box_per_row * box_per_row; // number of boxes on one face of the cube container
    
    // Clear out the box particle vector of old neighbors
    box_particles.clear();
    box_neighbors.clear();
    
    box_particles.reserve(box_per_row*box_per_row*box_per_row);
    box_neighbors.reserve(box_per_row*box_per_row*box_per_row);
    
    // Fill the box particles vector with empty vectors of ints
    for (int i = 0; i < box_per_row*box_per_row*box_per_row; i++) {
        box_particles.push_back(vector<int>());
        box_neighbors.push_back(vector<int>());
    }
    
    int box_num;
    for (int i = 0; i < num_particles; i++) {
        // determine box #
        box_num = compute_box_num(particles[i]->position, support_rad, max, min);
        // set box # in particle
        particles[i]->box = box_num;
        // add particle number to corresponding box
        box_particles[box_num].push_back(i);
    }
    
    // get all particles from neighboring boxes and add to particle's neighbor vector
    for (int i = 0; i < num_particles; i++) {
        particles[i]->neighbors.clear(); // clear out old neighbor vector
        int particle_num;
        vector<int> neighbor_boxes; // contains the numbers of all the neighboring boxes
        box_num = particles[i]->box;
        
        neighbor_boxes = surrounding_boxes(box_num, box_per_row, square_face);
        
        // add current particle as neighbor to current box
        int current_box_num = particles[i]->box;
        add_to_box_neighbors(current_box_num,i);
        
        // add particles in neighboring boxes to the current particles vector of neighbor particles
        // cycle through all neighboring boxes
        for (int j = 0; j < neighbor_boxes.size(); j++) {
            int num = neighbor_boxes[j];
            if (num < 0 || num > box_per_row*box_per_row*box_per_row) {
                cout<<"Error: incorrect neighboring box computed."<<endl;
                num = 0; // to prevent bad vector access
            }
            
            if (SURFACE) {
                add_to_box_neighbors(num,i);
            }
            
            // cycle through each particle in a box
            for (std::vector<int>::iterator it = box_particles[num].begin(); it != box_particles[num].end(); ++it) {
                particle_num = *it;
                Vec3 a = particles[i]->position;
                Vec3 b = particles[particle_num]->position;
                Vec3 diff = b-a;
                float dot_prod = dot(diff,diff);
                float dist = sqrt(dot_prod);
                
                if (SURFACE) {
                    add_to_box_neighbors(current_box_num,particle_num);
                }
                
                // check that the neighboring particle is within the support radius
                if (dist <= support_rad/2) {
                    particles[i]->neighbors.push_back(particle_num);
                }
            }
        }
    }
}

vector<int> Neighbor::surrounding_boxes(int box_num, int box_per_row, int square_face) {
    vector<int> neighbor_boxes;
    
    // each particle is a neighbor to particles in the same box
    neighbor_boxes.push_back(box_num);
    
    // add the correct numbers of neighboring boxes of the current box
    // 27 different cases
    if (box_num < square_face) {
        // back face
        if (box_num < box_per_row) {
            // bottom row
            if (box_num == 0) {
                // left side - back, bottom
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row+1);
                neighbor_boxes.push_back(box_num+square_face);
                neighbor_boxes.push_back(box_num+square_face+1);
                neighbor_boxes.push_back(box_num+square_face+box_per_row);
                neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
            } else if (box_num == box_per_row - 1) {
                // right side - back, bottom
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row-1);
                neighbor_boxes.push_back(box_num+square_face);
                neighbor_boxes.push_back(box_num+square_face-1);
                neighbor_boxes.push_back(box_num+square_face+box_per_row);
                neighbor_boxes.push_back(box_num+square_face+box_per_row-1);
            } else {
                // bottom, not left or right side - back, bottom
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row-1);
                neighbor_boxes.push_back(box_num+box_per_row+1);
                neighbor_boxes.push_back(box_num+square_face);
                neighbor_boxes.push_back(box_num+square_face-1);
                neighbor_boxes.push_back(box_num+square_face+1);
                neighbor_boxes.push_back(box_num+square_face+box_per_row);
                neighbor_boxes.push_back(box_num+square_face+box_per_row-1);
                neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
            }
        } else if (box_num >= square_face - box_per_row) {
            // top row
            if (box_num % box_per_row == 0) {
                // left side - back, top
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row+1);
                neighbor_boxes.push_back(box_num+square_face);
                neighbor_boxes.push_back(box_num+square_face+1);
                neighbor_boxes.push_back(box_num+square_face-box_per_row);
                neighbor_boxes.push_back(box_num+square_face-box_per_row+1);
                
            } else if (box_num % box_per_row == box_per_row - 1) {
                // right side - back, top
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row-1);
                neighbor_boxes.push_back(box_num+square_face);
                neighbor_boxes.push_back(box_num+square_face-1);
                neighbor_boxes.push_back(box_num+square_face-box_per_row);
                neighbor_boxes.push_back(box_num+square_face-box_per_row-1);
            } else {
                // top, not left or right - back, top
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row-1);
                neighbor_boxes.push_back(box_num-box_per_row+1);
                neighbor_boxes.push_back(box_num+square_face);
                neighbor_boxes.push_back(box_num+square_face-1);
                neighbor_boxes.push_back(box_num+square_face+1);
                neighbor_boxes.push_back(box_num+square_face-box_per_row);
                neighbor_boxes.push_back(box_num+square_face-box_per_row-1);
                neighbor_boxes.push_back(box_num+square_face-box_per_row+1);
            }
        } else if (box_num % box_per_row == 0) {
            // left side - back, not top or bottom
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row+1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row+1);
            neighbor_boxes.push_back(box_num+square_face);
            neighbor_boxes.push_back(box_num+square_face+1);
            neighbor_boxes.push_back(box_num+square_face+box_per_row);
            neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
            neighbor_boxes.push_back(box_num+square_face-box_per_row);
            neighbor_boxes.push_back(box_num+square_face-box_per_row+1);
            
        } else if (box_num % box_per_row == box_per_row - 1) {
            // right side - back, not top or bottom
            neighbor_boxes.push_back(box_num-1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row-1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row-1);
            neighbor_boxes.push_back(box_num+square_face);
            neighbor_boxes.push_back(box_num+square_face-1);
            neighbor_boxes.push_back(box_num+square_face+box_per_row);
            neighbor_boxes.push_back(box_num+square_face+box_per_row-1);
            neighbor_boxes.push_back(box_num+square_face-box_per_row);
            neighbor_boxes.push_back(box_num+square_face-box_per_row-1);
        } else {
            // back, not on edge
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num-1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row+1);
            neighbor_boxes.push_back(box_num+box_per_row-1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row+1);
            neighbor_boxes.push_back(box_num-box_per_row-1);
            neighbor_boxes.push_back(box_num+square_face);
            neighbor_boxes.push_back(box_num+square_face+1);
            neighbor_boxes.push_back(box_num+square_face-1);
            neighbor_boxes.push_back(box_num+square_face+box_per_row);
            neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
            neighbor_boxes.push_back(box_num+square_face+box_per_row-1);
            neighbor_boxes.push_back(box_num+square_face-box_per_row);
            neighbor_boxes.push_back(box_num+square_face-box_per_row+1);
            neighbor_boxes.push_back(box_num+square_face-box_per_row-1);
        }
        
    } else if (box_num >= square_face * (box_per_row - 1)) {
        // front face
        if (box_num % square_face < box_per_row) {
            // bottom row
            if (box_num % box_per_row == 0) {
                // left side - front, bottom
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row+1);
                neighbor_boxes.push_back(box_num-square_face);
                neighbor_boxes.push_back(box_num-square_face+1);
                neighbor_boxes.push_back(box_num-square_face+box_per_row);
                neighbor_boxes.push_back(box_num-square_face+box_per_row+1);
            } else if (box_num % box_per_row == box_per_row - 1) {
                // right side - front, bottom
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row-1);
                neighbor_boxes.push_back(box_num-square_face);
                neighbor_boxes.push_back(box_num-square_face-1);
                neighbor_boxes.push_back(box_num-square_face+box_per_row);
                neighbor_boxes.push_back(box_num-square_face+box_per_row-1);
            } else {
                // bottom, not left or right side - front, bottom
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num+box_per_row);
                neighbor_boxes.push_back(box_num+box_per_row+1);
                neighbor_boxes.push_back(box_num+box_per_row-1);
                neighbor_boxes.push_back(box_num-square_face);
                neighbor_boxes.push_back(box_num-square_face+1);
                neighbor_boxes.push_back(box_num-square_face-1);
                neighbor_boxes.push_back(box_num-square_face+box_per_row);
                neighbor_boxes.push_back(box_num-square_face+box_per_row+1);
                neighbor_boxes.push_back(box_num-square_face+box_per_row-1);
            }
        } else if (box_num % square_face >= square_face - box_per_row) {
            // top row
            if (box_num % box_per_row == 0) {
                // left side - front, top
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row+1);
                neighbor_boxes.push_back(box_num-square_face);
                neighbor_boxes.push_back(box_num-square_face+1);
                neighbor_boxes.push_back(box_num-square_face-box_per_row);
                neighbor_boxes.push_back(box_num-square_face-box_per_row+1);
            } else if (box_num % box_per_row == box_per_row - 1) {
                // right side - front, top
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row-1);
                neighbor_boxes.push_back(box_num-square_face);
                neighbor_boxes.push_back(box_num-square_face-1);
                neighbor_boxes.push_back(box_num-square_face-box_per_row);
                neighbor_boxes.push_back(box_num-square_face-box_per_row-1);
            } else {
                // top, not left or right - front, top
                neighbor_boxes.push_back(box_num+1);
                neighbor_boxes.push_back(box_num-1);
                neighbor_boxes.push_back(box_num-box_per_row);
                neighbor_boxes.push_back(box_num-box_per_row+1);
                neighbor_boxes.push_back(box_num-box_per_row-1);
                neighbor_boxes.push_back(box_num-square_face);
                neighbor_boxes.push_back(box_num-square_face+1);
                neighbor_boxes.push_back(box_num-square_face-1);
                neighbor_boxes.push_back(box_num-square_face-box_per_row);
                neighbor_boxes.push_back(box_num-square_face-box_per_row+1);
                neighbor_boxes.push_back(box_num-square_face-box_per_row-1);
            }
        } else if (box_num % box_per_row == 0) {
            // left side - front, not top or bottom
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row+1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row+1);
            neighbor_boxes.push_back(box_num-square_face);
            neighbor_boxes.push_back(box_num-square_face+1);
            neighbor_boxes.push_back(box_num-square_face+box_per_row);
            neighbor_boxes.push_back(box_num-square_face+box_per_row+1);
            neighbor_boxes.push_back(box_num-square_face-box_per_row);
            neighbor_boxes.push_back(box_num-square_face-box_per_row+1);
        } else if (box_num % box_per_row == box_per_row - 1) {
            // right side - front, not top or bottom
            neighbor_boxes.push_back(box_num-1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row-1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row-1);
            neighbor_boxes.push_back(box_num-square_face);
            neighbor_boxes.push_back(box_num-square_face-1);
            neighbor_boxes.push_back(box_num-square_face+box_per_row);
            neighbor_boxes.push_back(box_num-square_face+box_per_row-1);
            neighbor_boxes.push_back(box_num-square_face-box_per_row);
            neighbor_boxes.push_back(box_num-square_face-box_per_row-1);
        } else {
            // front face, not on edge
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num-1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row+1);
            neighbor_boxes.push_back(box_num+box_per_row-1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row+1);
            neighbor_boxes.push_back(box_num-box_per_row-1);
            neighbor_boxes.push_back(box_num-square_face);
            neighbor_boxes.push_back(box_num-square_face+1);
            neighbor_boxes.push_back(box_num-square_face-1);
            neighbor_boxes.push_back(box_num-square_face+box_per_row);
            neighbor_boxes.push_back(box_num-square_face+box_per_row+1);
            neighbor_boxes.push_back(box_num-square_face+box_per_row-1);
            neighbor_boxes.push_back(box_num-square_face-box_per_row);
            neighbor_boxes.push_back(box_num-square_face-box_per_row+1);
            neighbor_boxes.push_back(box_num-square_face-box_per_row-1);
        }
        
    } else if (box_num % box_per_row == 0) {
        // left face
        if (box_num % square_face == square_face - box_per_row) {
            // top - not front or back, left face
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row+1);
            neighbor_boxes.push_back(box_num-square_face);
            neighbor_boxes.push_back(box_num-square_face+1);
            neighbor_boxes.push_back(box_num-square_face-box_per_row);
            neighbor_boxes.push_back(box_num-square_face-box_per_row+1);
            neighbor_boxes.push_back(box_num+square_face);
            neighbor_boxes.push_back(box_num+square_face+1);
            neighbor_boxes.push_back(box_num+square_face-box_per_row);
            neighbor_boxes.push_back(box_num+square_face-box_per_row+1);
        } else if (box_num % square_face < box_per_row) {
            // bottom - not front or back, left face
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row+1);
            neighbor_boxes.push_back(box_num-square_face);
            neighbor_boxes.push_back(box_num-square_face+1);
            neighbor_boxes.push_back(box_num-square_face+box_per_row);
            neighbor_boxes.push_back(box_num-square_face+box_per_row+1);
            neighbor_boxes.push_back(box_num+square_face);
            neighbor_boxes.push_back(box_num+square_face+1);
            neighbor_boxes.push_back(box_num+square_face+box_per_row);
            neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
        } else {
            // left face - not on edge
            neighbor_boxes.push_back(box_num+1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row+1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row+1);
            neighbor_boxes.push_back(box_num-square_face);
            neighbor_boxes.push_back(box_num-square_face+1);
            neighbor_boxes.push_back(box_num-square_face+box_per_row);
            neighbor_boxes.push_back(box_num-square_face+box_per_row+1);
            neighbor_boxes.push_back(box_num-square_face-box_per_row);
            neighbor_boxes.push_back(box_num-square_face-box_per_row+1);
            neighbor_boxes.push_back(box_num+square_face);
            neighbor_boxes.push_back(box_num+square_face+1);
            neighbor_boxes.push_back(box_num+square_face+box_per_row);
            neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
            neighbor_boxes.push_back(box_num+square_face-box_per_row);
            neighbor_boxes.push_back(box_num+square_face-box_per_row+1);
        }
        
    } else if (box_num % box_per_row == box_per_row - 1) {
        // right face
        if (box_num % square_face == square_face - 1) {
            // top - not front or back, right face
            neighbor_boxes.push_back(box_num-1);
            neighbor_boxes.push_back(box_num-box_per_row);
            neighbor_boxes.push_back(box_num-box_per_row-1);
            neighbor_boxes.push_back(box_num-square_face);
            neighbor_boxes.push_back(box_num-square_face-1);
            neighbor_boxes.push_back(box_num-square_face-box_per_row);
            neighbor_boxes.push_back(box_num-square_face-box_per_row-1);
            neighbor_boxes.push_back(box_num+square_face);
            neighbor_boxes.push_back(box_num+square_face-1);
            neighbor_boxes.push_back(box_num+square_face-box_per_row);
            neighbor_boxes.push_back(box_num+square_face-box_per_row-1);
        } else if (box_num % square_face < box_per_row) {
            // bottom - not front or back, right face
            neighbor_boxes.push_back(box_num-1);
            neighbor_boxes.push_back(box_num+box_per_row);
            neighbor_boxes.push_back(box_num+box_per_row-1);
            neighbor_boxes.push_back(box_num-square_face);
            neighbor_boxes.push_back(box_num-square_face-1);
            neighbor_boxes.push_back(box_num-square_face+box_per_row);
            neighbor_boxes.push_back(box_num-square_face+box_per_row-1);
            neighbor_boxes.push_back(box_num+square_face);
            neighbor_boxes.push_back(box_num+square_face-1);
            neighbor_boxes.push_back(box_num+square_face+box_per_row);
            neighbor_boxes.push_back(box_num+square_face+box_per_row-1);
            neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
        }
        
        
        // add particles in neighboring boxes to the current particles vector of neighbor particles
        for (int j = 0; j < neighbor_boxes.size(); j++) {
            int num = neighbor_boxes[j];
            if (num < 0 || num > box_per_row*box_per_row*box_per_row) {
                cout<<"Error: incorrect neighboring box computed."<<endl;
                num = 0; // to prevent bad vector access
            }
            for (std::vector<int>::iterator it = box_particles[num].begin(); it != box_particles[num].end(); ++it) {
                particle_num = *it;
                Vec3 a = particles[i]->position;
                Vec3 b = particles[particle_num]->position;
				Vec3 diff = b-a;
                float dist = sqrt(dot(diff,diff));
                // check that the neighboring particle is within the support radius
                if (dist <= support_rad) {
                    particles[i]->neighbors.push_back(particle_num);
                }
            }
        }
    } else if (box_num % square_face >= square_face - box_per_row) {
        // top face
        neighbor_boxes.push_back(box_num+1);
        neighbor_boxes.push_back(box_num-1);
        neighbor_boxes.push_back(box_num-box_per_row);
        neighbor_boxes.push_back(box_num-box_per_row-1);
        neighbor_boxes.push_back(box_num-box_per_row+1);
        neighbor_boxes.push_back(box_num-square_face);
        neighbor_boxes.push_back(box_num-square_face-1);
        neighbor_boxes.push_back(box_num-square_face+1);
        neighbor_boxes.push_back(box_num-square_face-box_per_row);
        neighbor_boxes.push_back(box_num-square_face-box_per_row-1);
        neighbor_boxes.push_back(box_num-square_face-box_per_row+1);
        neighbor_boxes.push_back(box_num+square_face);
        neighbor_boxes.push_back(box_num+square_face-1);
        neighbor_boxes.push_back(box_num+square_face+1);
        neighbor_boxes.push_back(box_num+square_face-box_per_row);
        neighbor_boxes.push_back(box_num+square_face-box_per_row-1);
        neighbor_boxes.push_back(box_num+square_face-box_per_row+1);
    } else if (box_num % square_face < box_per_row) {
        // bottom face
        neighbor_boxes.push_back(box_num+1);
        neighbor_boxes.push_back(box_num-1);
        neighbor_boxes.push_back(box_num+box_per_row);
        neighbor_boxes.push_back(box_num+box_per_row-1);
        neighbor_boxes.push_back(box_num+box_per_row+1);
        neighbor_boxes.push_back(box_num-square_face);
        neighbor_boxes.push_back(box_num-square_face-1);
        neighbor_boxes.push_back(box_num-square_face+1);
        neighbor_boxes.push_back(box_num-square_face+box_per_row);
        neighbor_boxes.push_back(box_num-square_face+box_per_row-1);
        neighbor_boxes.push_back(box_num-square_face+box_per_row+1);
        neighbor_boxes.push_back(box_num+square_face);
        neighbor_boxes.push_back(box_num+square_face-1);
        neighbor_boxes.push_back(box_num+square_face+1);
        neighbor_boxes.push_back(box_num+square_face+box_per_row);
        neighbor_boxes.push_back(box_num+square_face+box_per_row-1);
        neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
    } else {
        // all other particles, inside cube
        neighbor_boxes.push_back(box_num+1);
        neighbor_boxes.push_back(box_num-1);
        neighbor_boxes.push_back(box_num-box_per_row);
        neighbor_boxes.push_back(box_num-box_per_row-1);
        neighbor_boxes.push_back(box_num-box_per_row+1);
        neighbor_boxes.push_back(box_num+box_per_row);
        neighbor_boxes.push_back(box_num+box_per_row-1);
        neighbor_boxes.push_back(box_num+box_per_row+1);
        neighbor_boxes.push_back(box_num-square_face);
        neighbor_boxes.push_back(box_num-square_face-1);
        neighbor_boxes.push_back(box_num-square_face+1);
        neighbor_boxes.push_back(box_num-square_face-box_per_row);
        neighbor_boxes.push_back(box_num-square_face-box_per_row-1);
        neighbor_boxes.push_back(box_num-square_face-box_per_row+1);
        neighbor_boxes.push_back(box_num-square_face+box_per_row);
        neighbor_boxes.push_back(box_num-square_face+box_per_row-1);
        neighbor_boxes.push_back(box_num-square_face+box_per_row+1);
        neighbor_boxes.push_back(box_num+square_face);
        neighbor_boxes.push_back(box_num+square_face-1);
        neighbor_boxes.push_back(box_num+square_face+1);
        neighbor_boxes.push_back(box_num+square_face-box_per_row);
        neighbor_boxes.push_back(box_num+square_face-box_per_row-1);
        neighbor_boxes.push_back(box_num+square_face-box_per_row+1);
        neighbor_boxes.push_back(box_num+square_face+box_per_row);
        neighbor_boxes.push_back(box_num+square_face+box_per_row-1);
        neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
    }
    return neighbor_boxes;
}