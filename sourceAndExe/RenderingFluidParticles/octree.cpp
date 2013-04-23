
#include "FluidSimulation.h"

#define MINIMUM_DEPTH 6
#define MAX_DEPTH 9
#define EXTRACT_ISOVALUE 0.01



void Octree::InitOctree(vector<WaterParticle> &water_particles)
{
    for(int i = 0; i < 3; i++) root.bbmin[i] =  1.0e6;
    for(int i = 0; i < 3; i++) root.bbmax[i] = -1.0e6;

    for(unsigned int i = 0; i < water_particles.size(); i++){
        for(int j = 0; j < 3; j++){
            if(water_particles[i].position[j] < root.bbmin[j]) root.bbmin[j] = water_particles[i].position[j];
            if(water_particles[i].position[j] > root.bbmax[j]) root.bbmax[j] = water_particles[i].position[j];
        }
    }

    int largest_dim = -1;

    if(root.bbmax[0]-root.bbmin[0] > root.bbmax[1]-root.bbmin[1]){
        if(root.bbmax[0]-root.bbmin[0] > root.bbmax[2]-root.bbmin[2]) 
            largest_dim = 0; 
        else
            largest_dim = 2;
    }else{
        if(root.bbmax[1]-root.bbmin[1] > root.bbmax[2]-root.bbmin[2]) 
            largest_dim = 1; 
        else
            largest_dim = 2;
    }

    double largest_dim_size = (root.bbmax[largest_dim]-root.bbmin[largest_dim]) + H*2.0;

    double center[3];
    for(int j = 0; j < 3; j++) center[j] = (root.bbmin[j]+root.bbmax[j])*0.5;
    
    for(int j = 0; j < 3; j++){
        root.bbmin[j] = center[j] - largest_dim_size*0.5;
        root.bbmax[j] = center[j] + largest_dim_size*0.5;
    }

    ConstructOctree(&root, water_particles, 0);
}



void Octree::ConstructOctree(OctreeNode *node, vector<WaterParticle> &water_particles, int depth)
{

    if(depth >= MINIMUM_DEPTH){

        for(int i = 0; i < 3; i++) node->inv_size[i] = 1.0 / (node->bbmax[i]-node->bbmin[i]);

        node->surface_distance.resize(8);

        for(int id = 0; id < 8; id++){

            double cell_vertex[3];

            switch(id){
            case 0: 
                cell_vertex[0] = node->bbmin[0];  cell_vertex[1] = node->bbmin[1];  cell_vertex[2] = node->bbmin[2];
                break;
            case 1: 
                cell_vertex[0] = node->bbmax[0];  cell_vertex[1] = node->bbmin[1];  cell_vertex[2] = node->bbmin[2];
                break;
            case 2: 
                cell_vertex[0] = node->bbmax[0];  cell_vertex[1] = node->bbmax[1];  cell_vertex[2] = node->bbmin[2];
                break;
            case 3: 
                cell_vertex[0] = node->bbmin[0];  cell_vertex[1] = node->bbmax[1];  cell_vertex[2] = node->bbmin[2];
                break;
            case 4: 
                cell_vertex[0] = node->bbmin[0];  cell_vertex[1] = node->bbmin[1];  cell_vertex[2] = node->bbmax[2];
                break;
            case 5: 
                cell_vertex[0] = node->bbmax[0];  cell_vertex[1] = node->bbmin[1];  cell_vertex[2] = node->bbmax[2];
                break;
            case 6: 
                cell_vertex[0] = node->bbmax[0];  cell_vertex[1] = node->bbmax[1];  cell_vertex[2] = node->bbmax[2];
                break;
            case 7: 
                cell_vertex[0] = node->bbmin[0];  cell_vertex[1] = node->bbmax[1];  cell_vertex[2] = node->bbmax[2];
                break;
            }

            double center_of_mass[3] = { 0.0, 0.0, 0.0, };
            double sumW = 0.0;

            for(unsigned int i = 0; i < water_particles.size(); i++){

                double r[3] = { cell_vertex[0] - water_particles[i].position[0], 
                                cell_vertex[1] - water_particles[i].position[1], 
                                cell_vertex[2] - water_particles[i].position[2] };
                double r_square = DotProduct(r, r);

                static double support_r_square = H*H;

                if(r_square <= support_r_square){                    
                    double w_poly6 = Wpoly6(r_square);

                    for(int j = 0; j < 3; j++) center_of_mass[j] += water_particles[i].position[j]*w_poly6;
                    sumW += w_poly6;
                } 

            }

            double difference_vec[3];
            double inv_sumW = 1.0/sumW;

            for(int j = 0; j < 3; j++) difference_vec[j] = cell_vertex[j] - center_of_mass[j]*inv_sumW;

            node->surface_distance[id] = sqrt( DotProduct(difference_vec, difference_vec) );

            // no neighbor particles
            if(fabs(sumW) < 1.0e-12) node->surface_distance[id] = 1.0e6;

        } // for(int id = 0; id < 8; id++){

        bool is_below_extract_isovalue = false, is_above_extract_isovalue = false;


        for(int id = 0; id < 8; id++){
            if(node->surface_distance[id] < EXTRACT_ISOVALUE) is_below_extract_isovalue = true;
            if(node->surface_distance[id] > EXTRACT_ISOVALUE) is_above_extract_isovalue = true;
        }

        if(depth == MAX_DEPTH || !is_below_extract_isovalue || !is_above_extract_isovalue ) 
            return;

     //   return;
    }


 
    for(int id = 0; id < 8; id++){
        if(node->child[id] == NULL) AddChild(node, id);

        GetChildBBRange(node->bbmin, node->bbmax, node->child[id]->bbmin, node->child[id]->bbmax, id);
    }


    vector<WaterParticle> water_particles_child_node[8];

    for(unsigned int i = 0; i < water_particles.size(); i++){

        double *pos = water_particles[i].position;

        for(int id = 0; id < 8; id++){
            if( node->child[id]->bbmin[0]-H < pos[0] && pos[0] < node->child[id]->bbmax[0]+H && 
                node->child[id]->bbmin[1]-H < pos[1] && pos[1] < node->child[id]->bbmax[1]+H &&
                node->child[id]->bbmin[2]-H < pos[2] && pos[2] < node->child[id]->bbmax[2]+H ){

                water_particles_child_node[id].push_back(water_particles[i]);
            }
        }

    } // vector<WaterParticle> water_particles_child_node[8]


    for(int id = 0; id < 8; id++){
        if(water_particles_child_node[id].empty() == false){
            ConstructOctree(node->child[id], water_particles_child_node[id], depth+1);
        }
    }

}

void Octree::AddChild(OctreeNode *node, int id)
{
  OctreeNode *new_node;

  new_node = new(nothrow) OctreeNode;
  if(!new_node){ cerr << "memory allocation failed\n"; exit(-1); }

  node->child[id]  = new_node;
  new_node->parent = node;

  for(int i = 0; i < 8; i++) new_node->child[i] = NULL;
}

void Octree::DeleteOctree()
{
    ReleaseNodeMemory(&root);
}


void Octree::ReleaseNodeMemory(OctreeNode *node)
{
    // leaf
    if(node->child[0] == NULL) return;
    
    for(int id = 0; id < 8; id++) ReleaseNodeMemory(node->child[id]);
    
    for(int i = 0; i < 8; i++){
        delete (node->child[i]);
        node->child[i] = NULL;
    }
}


void Octree::GetChildBBRange(double *bbmin, double *bbmax, double *new_bbmin, double *new_bbmax, int id)
{
    double mid[3];

    for(int i = 0; i < 3; i++)
        mid[i] = (bbmin[i] + bbmax[i]) * 0.5;

    if(id < 4){

        if(id < 2){

            if(id == 0){
                new_bbmin[0] = bbmin[0];
                new_bbmax[0] = mid[0];
            }else{ // id == 1
                new_bbmin[0] = mid[0];
                new_bbmax[0] = bbmax[0];
            }

            new_bbmin[1] = bbmin[1];
            new_bbmax[1] = mid[1];

        }else{

            if(id == 2){
                new_bbmin[0] = bbmin[0];
                new_bbmax[0] = mid[0];
            }else{ // id == 3
                new_bbmin[0] = mid[0];
                new_bbmax[0] = bbmax[0];
            }

            new_bbmin[1] = mid[1];
            new_bbmax[1] = bbmax[1];

        }
    
        new_bbmin[2] = bbmin[2];
        new_bbmax[2] = mid[2];
   
    }else{

        if(id < 6){

            if(id == 4){
                new_bbmin[0] = bbmin[0];
                new_bbmax[0] = mid[0];
            }else{ // id == 5
                new_bbmin[0] = mid[0];
                new_bbmax[0] = bbmax[0];
            }

            new_bbmin[1] = bbmin[1];
            new_bbmax[1] = mid[1];

        }else{

            if(id == 6){
                new_bbmin[0] = bbmin[0];
                new_bbmax[0] = mid[0];
            }else{ // id == 7
                new_bbmin[0] = mid[0];
                new_bbmax[0] = bbmax[0];
            }

            new_bbmin[1] = mid[1];
            new_bbmax[1] = bbmax[1];
        }

        new_bbmin[2] = mid[2];
        new_bbmax[2] = bbmax[2];
    }

}


double TrilinearInterpolation(OctreeNode *node, double u, double v, double w)
{
    double a = node->surface_distance[0];
    double b = node->surface_distance[1];
    double c = node->surface_distance[2];
    double d = node->surface_distance[3];
    double e = node->surface_distance[4];
    double f = node->surface_distance[5];
    double g = node->surface_distance[6];
    double h = node->surface_distance[7];

    return ( (h*w + d*(1-w))*v + (e*w + a*(1-w))*(1-v) )*(1-u) + ( (g*w + c*(1-w))*v + (f*w + b*(1-w))*(1-v) )*u;
}


void GetTrilinearParameter(OctreeNode *node, double *point, double &u, double &v, double &w)
{
    u = (point[0] - node->bbmin[0]) * node->inv_size[0];
    v = (point[1] - node->bbmin[1]) * node->inv_size[1];
    w = (point[2] - node->bbmin[2]) * node->inv_size[2];
}

bool IsRayIntersectWaterSurface(OctreeNode *node, double *entry, double *exit, double *intersection, double *gradient, bool isComputeIntersection)
{
    double u_in,  v_in,  w_in; 
    double u_out, v_out, w_out;

    GetTrilinearParameter(node, entry,  u_in,  v_in,  w_in);
    GetTrilinearParameter(node, exit,  u_out, v_out, w_out);

    double color_field_in  = TrilinearInterpolation(node,  u_in,  v_in,  w_in);
    double color_field_out = TrilinearInterpolation(node, u_out, v_out, w_out);


    if( (color_field_in < EXTRACT_ISOVALUE && color_field_out > EXTRACT_ISOVALUE) ||  
        (color_field_in > EXTRACT_ISOVALUE && color_field_out < EXTRACT_ISOVALUE) ){

        if(isComputeIntersection){
            double intersect_t = (EXTRACT_ISOVALUE - color_field_in) / (color_field_out - color_field_in);

            double u_hit = u_in*(1.0-intersect_t) + u_out*intersect_t;
            double v_hit = v_in*(1.0-intersect_t) + v_out*intersect_t;
            double w_hit = w_in*(1.0-intersect_t) + w_out*intersect_t;

            gradient[0] = TrilinearInterpolation(node, 1.0, v_hit, w_hit) - TrilinearInterpolation(node, 0.0, v_hit, w_hit);
            gradient[1] = TrilinearInterpolation(node, u_hit, 1.0, w_hit) - TrilinearInterpolation(node, u_hit, 0.0, w_hit);
            gradient[2] = TrilinearInterpolation(node, u_hit, v_hit, 1.0) - TrilinearInterpolation(node, u_hit, v_hit, 0.0);
            
            Normalize(gradient);

            for(int i = 0; i < 3; i++) intersection[i] = entry[i]*(1.0-intersect_t) + exit[i]*intersect_t;
        }

        return true;
    }else{
        return false;
    }
}


void Octree::ClipRaySegment(double *ray_origin, double *ray_direction, double *ray_direction_inv, double *bbmin, double *bbmax,
                            double &t_near, double &t_far)
{
    t_near = EPSILON;
    t_far  = 1.0e6;

    for(int i = 0; i < 3; i++){

        if( fabs(ray_direction_inv[i]) > EPSILON ){

            double t1 = (bbmin[i]-ray_origin[i]) * ray_direction_inv[i];
            double t2 = (bbmax[i]-ray_origin[i]) * ray_direction_inv[i];      

            if(ray_direction_inv[i] > 0.0){
                if(t1 > t_near) t_near = t1;
                if(t2 < t_far)  t_far  = t2;
            }else{
                if(t1 < t_far)  t_far  = t1;        
                if(t2 > t_near) t_near = t2;
            }

        } // if( fabs(ray_direction_inv[i]) > EPSILON ){

    } // for(int i = 0; i < 3; i++){

}



bool Octree::FindIntersection(double *org, double *dir, double *inv_dir, double *intersection, double *gradient, bool isComputeIntersection)
{
    double t_near, t_far;
    ClipRaySegment(org, dir, inv_dir, root.bbmin, root.bbmax, t_near, t_far);


    if(t_near > t_far) return false; // ray does not hit root

    return CheckNode(&root, org, dir, inv_dir, t_near, t_far, intersection, gradient, isComputeIntersection);
}



bool Octree::CheckNode(OctreeNode *node, double *org, double *dir, double *inv_dir, double t_near, double t_far, 
                       double *intersection, double *gradient, bool isComputeIntersection)
{
    if(node->child[0] == NULL){
        if(node->surface_distance.size() == 0) return false;

        double entry[3], exit[3];

        for(int i = 0; i < 3; i++) entry[i] = org[i] + dir[i]*t_near;
        for(int i = 0; i < 3; i++) exit[i]  = org[i] + dir[i]*t_far;

        return IsRayIntersectWaterSurface(node, entry, exit, intersection, gradient, isComputeIntersection);
    }

    double center[3];
    for(int j = 0; j < 3; j++) center[j] = (node->bbmin[j]+node->bbmax[j])*0.5;


    double hit_point[3];
    for(int i = 0; i < 3; i++) hit_point[i] = org[i]+dir[i]*(t_near + 1.0e-6);

    while(node->bbmin[0] < hit_point[0] && hit_point[0] < node->bbmax[0] && 
          node->bbmin[1] < hit_point[1] && hit_point[1] < node->bbmax[1] && 
          node->bbmin[2] < hit_point[2] && hit_point[2] < node->bbmax[2]  ){

        int id = 0;

        if(hit_point[0] > center[0]) id++;
        if(hit_point[1] > center[1]) id += 2;
        if(hit_point[2] > center[2]) id += 4;

        double t_near_child, t_far_child;
        ClipRaySegment(org, dir, inv_dir, node->child[id]->bbmin, node->child[id]->bbmax, t_near_child, t_far_child);

        if( CheckNode(node->child[id], org, dir, inv_dir, t_near_child, t_far_child, intersection, gradient, isComputeIntersection) )
            return true;

        for(int i = 0; i < 3; i++) hit_point[i] = org[i]+dir[i]*(t_far_child + 1.0e-6);
    }

    return false;
}


