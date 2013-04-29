
#include "Neighbor.h"

#include <cmath>

#include <iostream>

#include <algorithm>

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


int Neighbor::compute_box_num(Vec3 pos, float support_rad, Vec3 max_point, Vec3 min_point) {
    // assuming container is rectangular, axis aligned
    int x = -1,y = -1, z = -1;
    float x_span = max_point.x - min_point.x;
    float y_span = max_point.y - min_point.y;
    float z_span = max_point.z - min_point.z;
    float max_span = max(x_span,y_span);
    max_span = max(max_span,z_span);
    int x_divs = (int)(x_span / support_rad);
    int y_divs = (int)(y_span / support_rad);
    int z_divs = (int)(z_span / support_rad);
    // int box_per_row = (int)max_span / support_rad; // casting to int, assuming support radius evenly divides width
    
    float curr_x = min_point.x, curr_y = min_point.y, curr_z = min_point.z;
    
	float col_point, row_point, depth_point;
    
	////if(col_point>curr_x || row_point>curr_y || col_point<0 || row_point<0){
	////	point is not inside the container.
	////	cout<<'h'<<endl;
	////	return 0;
	////}
    
	////row = floor(row_point/support_rad);
	////col = floor(col_point/support_rad);
    
	//int num = col + row*box_per_row - 1;
	
    ////return num;//int(max(float(num),0.0f));
    
	for (int i = 0; i < max_span && (x == -1 || y == -1 || z == -1); i++) {
        col_point = abs(pos.x - curr_x);
        row_point = abs(pos.y - curr_y);
        depth_point = abs(pos.z - curr_z);
        
		if (col_point <= support_rad && x == -1) {
			x = i;
		}
        if (row_point <= support_rad && y == -1) {
			y = i;
		}
        if (depth_point <= support_rad && z == -1) {
            z = i;
        }
        
        curr_x += support_rad;
        curr_y += support_rad;
        curr_z += support_rad;
	}
    
    int num = x + y * y_span + z * y_span * z_span;
    //cout<<"Num = "<<num<<" = "<<col<<" + "<<row<<" * "<<box_per_row<<endl;
    // cout<<endl;
    
    if (num >= x_span * y_span * z_span || num < 0 || x == -1 || y == -1 || z == -1) {
        //cout<<"Error, incorrect box # assigned in Neighbor: "<<num<<endl;
        num = 0; // set box number to 0 to prevent bad vector access
        // exit(0);
    }
    //compute box num given row & col
    return num;
}

void Neighbor::place_particles(vector<Particle*> &particles, float support_rad, Container c) {
    // assuming square container
//    float width = c.max.x - c.min.x;
    float x_span = c.max.x - c.min.x;
    float y_span = c.max.y - c.min.y;
    float z_span = c.max.z - c.min.z;
    
    int box_per_row = x_span/support_rad; // change
    int square_face = box_per_row * box_per_row;
    
    box_particles.clear();
    
    for (int i = 0; i < box_per_row*box_per_row*box_per_row; i++) {
        box_particles.push_back(vector<int>());
    }
    
    int box_num;
    for (int i = 0; i < particles.size(); i++) {
        // determine box #
        box_num = compute_box_num(particles[i]->position, support_rad, c.max, c.min);
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
        
        if (box_num < square_face) {
            // front face
            if (box_num < box_per_row) {
                // bottom row
                if (box_num == 0) {
                    // left side - front, bottom
                    neighbor_boxes.push_back(box_num+1);
                    neighbor_boxes.push_back(box_num+box_per_row);
                    neighbor_boxes.push_back(box_num+box_per_row+1);
                    neighbor_boxes.push_back(box_num+square_face);
                    neighbor_boxes.push_back(box_num+square_face+1);
                    neighbor_boxes.push_back(box_num+square_face+box_per_row);
                    neighbor_boxes.push_back(box_num+square_face+box_per_row+1);
                } else if (box_num == box_per_row - 1) {
                    // right side - front, bottom
                    neighbor_boxes.push_back(box_num-1);
                    neighbor_boxes.push_back(box_num+box_per_row);
                    neighbor_boxes.push_back(box_num+box_per_row-1);
                    neighbor_boxes.push_back(box_num+square_face);
                    neighbor_boxes.push_back(box_num+square_face-1);
                    neighbor_boxes.push_back(box_num+square_face+box_per_row);
                    neighbor_boxes.push_back(box_num+square_face+box_per_row-1);
                } else {
                    // bottom, not left or right side - front, bottom
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
                    // left side - front, top
                    neighbor_boxes.push_back(box_num+1);
                    neighbor_boxes.push_back(box_num-box_per_row);
                    neighbor_boxes.push_back(box_num-box_per_row+1);
                    neighbor_boxes.push_back(box_num+square_face);
                    neighbor_boxes.push_back(box_num+square_face+1);
                    neighbor_boxes.push_back(box_num+square_face-box_per_row);
                    neighbor_boxes.push_back(box_num+square_face-box_per_row+1);
                    
                } else if (box_num % box_per_row == box_per_row - 1) {
                    // right side - front, top
                    neighbor_boxes.push_back(box_num-1);
                    neighbor_boxes.push_back(box_num-box_per_row);
                    neighbor_boxes.push_back(box_num-box_per_row-1);
                    neighbor_boxes.push_back(box_num+square_face);
                    neighbor_boxes.push_back(box_num+square_face-1);
                    neighbor_boxes.push_back(box_num+square_face-box_per_row);
                    neighbor_boxes.push_back(box_num+square_face-box_per_row-1);
                } else {
                    // top, not left or right - front, top
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
                // left side - front, not top or bottom
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
                // right side - front, not top or bottom
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
                // front, not on edge
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
            // back face
            if (box_num % square_face < box_per_row) {
                // bottom row
                if (box_num % box_per_row == 0) {
                    // left side - back, bottom
                    neighbor_boxes.push_back(box_num+1);
                    neighbor_boxes.push_back(box_num+box_per_row);
                    neighbor_boxes.push_back(box_num+box_per_row+1);
                    neighbor_boxes.push_back(box_num-square_face);
                    neighbor_boxes.push_back(box_num-square_face+1);
                    neighbor_boxes.push_back(box_num-square_face+box_per_row);
                    neighbor_boxes.push_back(box_num-square_face+box_per_row+1);
                } else if (box_num % box_per_row == box_per_row - 1) {
                    // right side - back, bottom
                    neighbor_boxes.push_back(box_num-1);
                    neighbor_boxes.push_back(box_num+box_per_row);
                    neighbor_boxes.push_back(box_num+box_per_row-1);
                    neighbor_boxes.push_back(box_num-square_face);
                    neighbor_boxes.push_back(box_num-square_face-1);
                    neighbor_boxes.push_back(box_num-square_face+box_per_row);
                    neighbor_boxes.push_back(box_num-square_face+box_per_row-1);
                } else {
                    // bottom, not left or right side - back, bottom
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
                    // left side - back, top
                    neighbor_boxes.push_back(box_num+1);
                    neighbor_boxes.push_back(box_num-box_per_row);
                    neighbor_boxes.push_back(box_num-box_per_row+1);
                    neighbor_boxes.push_back(box_num-square_face);
                    neighbor_boxes.push_back(box_num-square_face+1);
                    neighbor_boxes.push_back(box_num-square_face-box_per_row);
                    neighbor_boxes.push_back(box_num-square_face-box_per_row+1);
                } else if (box_num % box_per_row == box_per_row - 1) {
                    // right side - back, top
                    neighbor_boxes.push_back(box_num-1);
                    neighbor_boxes.push_back(box_num-box_per_row);
                    neighbor_boxes.push_back(box_num-box_per_row-1);
                    neighbor_boxes.push_back(box_num-square_face);
                    neighbor_boxes.push_back(box_num-square_face-1);
                    neighbor_boxes.push_back(box_num-square_face-box_per_row);
                    neighbor_boxes.push_back(box_num-square_face-box_per_row-1);
                } else {
                    // top, not left or right - back, top
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
                // left side - back, not top or bottom
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
                // right side - back, not top or bottom
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
                // back face, not on edge
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
            } else if (box_num % box_per_row == 0) {
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
            } else if (box_num % box_per_row == box_per_row - 1) {
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
            } else {
                // right face - not on edge
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
                neighbor_boxes.push_back(box_num+square_face);
                neighbor_boxes.push_back(box_num+square_face-1);
                neighbor_boxes.push_back(box_num+square_face+box_per_row);
                neighbor_boxes.push_back(box_num+square_face+box_per_row-1);
                neighbor_boxes.push_back(box_num+square_face-box_per_row);
                neighbor_boxes.push_back(box_num+square_face-box_per_row-1);
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
        
        
        // need to find out which boxes neighbor the current box
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
                float dist = sqrt(pow((a.x - b.x),2) - pow((a.y - b.y),2));
                if (dist <= support_rad) {
                    particles[i]->neighbors.push_back(particle_num);
                }
            }
        }
    }
}
