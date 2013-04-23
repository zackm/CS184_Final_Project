
#include "FluidSimulation.h"


#define EXTRACT_ISOVALUE 0.01
#define MAX_N_ITERATION 5

#define AMBIENT 0.1
#define N_LIGHT 2

double light_pos[2][3] = { {-1, 0.7, 1},  { 1, 0.7, 1} };



int SPH::FindNearestObject(double *org, double *dir, double *inv_dir, double *hit_point, double *normal)
{
    int    nearest_object_index = -1;
    double nearest_t = 1.0e6+1;
   

    for(int i = 0; i < 5; i++){
        // if(i==3) continue;

        double bbmin[3], bbmax[3];

        switch(i){
        case 0: 
            bbmin[0] = -0.45;  bbmin[1] = 0.18;  bbmin[2] = 0.0;  
            bbmax[0] = -0.27;  bbmax[1] = 0.90;  bbmax[2] = 0.18;
            break;
        case 1: 
            bbmin[0] = -0.09;  bbmin[1] = 0.00;  bbmin[2] = 0.0;  
            bbmax[0] =  0.09;  bbmax[1] = 0.36;  bbmax[2] = 0.18;
            break;
        case 2: 
            bbmin[0] =  0.27;  bbmin[1] = 0.18;  bbmin[2] = 0.0;  
            bbmax[0] =  0.45;  bbmax[1] = 0.90;  bbmax[2] = 0.18;
            break;
        case 3: 
            bbmin[0] =  -0.27;  bbmin[1] = 0.72;  bbmin[2] = 0.0;  
            bbmax[0] =   0.27;  bbmax[1] = 0.90;  bbmax[2] = 0.18;
            break;
        case 4: // container glass surface
            bbmin[0] = -0.64;  bbmin[1] = -0.1;  bbmin[2] = -0.01;  
            bbmax[0] =  0.64;  bbmax[1] = 0.90;  bbmax[2] = 0.19;
            break;
        }

        double local_t;
        double local_intersection[3], local_normal[3];

        if( ComputeRayBoxIntersection(org, dir, inv_dir, bbmin, bbmax, local_t, local_intersection, local_normal) ){
            if(local_t < nearest_t){
                nearest_t = local_t;
                nearest_object_index = i;

                for(int j = 0; j < 3; j++) hit_point[j] = local_intersection[j];
                for(int j = 0; j < 3; j++) normal[j]    = local_normal[j];
            }

        }

    }


    // intersection test between ray and floor
    if(org[1] > -0.0001 && dir[1] < 0.0){     
        double t = (-0.0001-org[1])*inv_dir[1]; 
        if(t > 1.0e-6 && t < nearest_t){
            nearest_t = t;
            nearest_object_index = 5;

            normal[0] = 0.0;  
            normal[1] = 1.0;  
            normal[2] = 0.0; 
        }
    }

    // intersection test between ray and back wall
    if(dir[2] < 0.0){     
        double t = (-0.5-org[2])*inv_dir[2]; 
        if(t > 1.0e-6 && t < nearest_t){
            nearest_t = t;
            nearest_object_index = 6;

            normal[0] = 0.0;  
            normal[1] = 0.0;  
            normal[2] = 1.0; 
        }
    }

    for(int i = 0; i < 3; i++) hit_point[i] = org[i] + dir[i]*nearest_t;


    double water_surface_intersection[3], gradient[3];

    double water_bbmin[3] = { -0.63, 0.0, 0.0 }, water_bbmax[3] = { 0.63, 0.90, 0.18 };
    double dummy, dummy_array[3];

    // first perform bounding box test to avoid unnecessary octree traversal
    if( ComputeRayBoxIntersection(org, dir, inv_dir, water_bbmin, water_bbmax, dummy, dummy_array, dummy_array) ){

        if( octree.FindIntersection(org, dir, inv_dir, water_surface_intersection, gradient, true) ){
            for(int i = 0; i < 3; i++){           
                if(fabs(inv_dir[i]) > 1.0e-6){    
                  
                    double diff_t = (water_surface_intersection[i]-hit_point[i])*inv_dir[i];
                    if(diff_t < 0.0){
                        // the closest object is water surface
                        nearest_object_index = 100;
                        for(int i = 0; i < 3; i++) hit_point[i] = water_surface_intersection[i];
                        for(int i = 0; i < 3; i++) normal[i]    = gradient[i];
                    }
                    break;

                }
            }
        } // if( octree.FindIntersection(org, dir, inv_dir, water_surface_intersection, gradient, true) ){

    }


    return nearest_object_index;
}



void SPH::ComputeDirectLighting(double *point, double *normal, double *direction, int object_index, double depth, double *color)
{
    color[0] = color[1] = color[2] = 0.0;//= AMBIENT; // for ambient


    for(int i = 0; i < N_LIGHT; i++){

        double shadow_ray[3] = { light_pos[i][0]-point[0], light_pos[i][1]-point[1], light_pos[i][2]-point[2], };
        Normalize(shadow_ray);

        double inv_shadow_ray[3] = { 1.0/shadow_ray[0], 1.0/shadow_ray[1], 1.0/shadow_ray[2], };

        double attenuation_weight = 1.0;
        TraceShadowRay(point, shadow_ray, inv_shadow_ray, attenuation_weight); 


        double local_color[3], reflect_color[3];

        double reflection_vector[3];
        ComputeReflectionVector(normal, direction, reflection_vector);

        if( (object_index == 0 || object_index == 1 || object_index == 2 || object_index == 3 ) && depth < 1){

            double inv_reflection_vector[3];
            for(int j = 0; j < 3; j++) inv_reflection_vector[j] = 1.0 / reflection_vector[j];

            double reflect_obj_hit_point[3], reflect_obj_normal[3];

            int reflect_obj_index = FindNearestObject(point, reflection_vector, inv_reflection_vector, reflect_obj_hit_point, reflect_obj_normal);

            if( (reflect_obj_index != -1 ) ){
                RayCasting(reflect_obj_hit_point, reflection_vector,inv_reflection_vector, depth+1, reflect_color);
            }

        }else{
            reflect_color[0] = reflect_color[1] = reflect_color[2] = 0.0;
        }

        double diffuse_weight  = fabs( DotProduct(normal, shadow_ray) );
        double specular_weight = fabs( DotProduct(normal, reflection_vector) );


        if(object_index >= 0 && object_index <= 3){
            // container surface
            pow(specular_weight, 10);

            local_color[0] = diffuse_weight*0.7 + specular_weight*0.3*reflect_color[0];
            local_color[1] = diffuse_weight*0.7 + specular_weight*0.3*reflect_color[1];
            local_color[2] = diffuse_weight*0.7 + specular_weight*0.3*reflect_color[2];

        }else if(object_index == 4){
            // container surface
            specular_weight = pow(specular_weight, 10);

            local_color[0] = diffuse_weight*0.5 + specular_weight*0.5;//*reflect_color[0];
            local_color[1] = diffuse_weight*0.5 + specular_weight*0.5;//*reflect_color[1];
            local_color[2] = diffuse_weight*0.5 + specular_weight*0.5;//*reflect_color[2];

        }else if(object_index == 5){
            // floor

            if(point[0] > -1.0 && point[0] < 1.0 && point[2] > -0.5 && point[2] < 1.5 ){
                double tex_coord[2] = { (point[0]+1.0)*0.5, (point[2]+0.5)*0.5 };

                int pixel_index =  png_width*(int)(tex_coord[1]*png_height) + tex_coord[0]*png_width;               

                local_color[0] = diffuse_weight * png_data[pixel_index*3];
                local_color[1] = diffuse_weight * png_data[pixel_index*3+1];
                local_color[2] = diffuse_weight * png_data[pixel_index*3+2];

            }else{
                local_color[0] = local_color[1] = local_color[2] = diffuse_weight*1.0;
            }

        }else if(object_index == 6){
            // back wall
            specular_weight = pow(specular_weight, 10);

            local_color[0] = diffuse_weight*0.7 + specular_weight*0.0;
            local_color[1] = diffuse_weight*0.7 + specular_weight*0.0;
            local_color[2] = diffuse_weight*0.4 + specular_weight*0.0;

        }else if(object_index == 100){
            // water surface
            specular_weight = pow(specular_weight, 10);

            local_color[0] = diffuse_weight*0.0 + specular_weight*0.6;
            local_color[1] = diffuse_weight*0.1 + specular_weight*0.6;
            local_color[2] = diffuse_weight*0.1 + specular_weight*0.6;

        }

        color[0] += attenuation_weight * local_color[0];
        color[1] += attenuation_weight * local_color[1];
        color[2] += attenuation_weight * local_color[2];

    } // for(int i = 0; i < 2; i++){

}




void SPH::TraceShadowRay(double *org, double *dir, double *inv_dir, double &weight)
{
    static double attenuation_rate_container = 0.9;
    static double attenuation_rate_water     = 0.7;

    double hit_point[3], normal[3];

    int nearest_object_index = FindNearestObject(org, dir, inv_dir, hit_point, normal);
    

    if(nearest_object_index == -1){ // hit light

        return;
            
    }else if(nearest_object_index >=0 && nearest_object_index <= 3){
        
        weight *= 0.2;
        return;

    }else if(nearest_object_index == 4){
        
        weight *= attenuation_rate_container; 
        TraceShadowRay(hit_point, dir, inv_dir, weight);
        return;

    }else if(nearest_object_index == 5 || nearest_object_index == 6){
        
//        weight = 1.0;
        return;

    }else if(nearest_object_index == 100){ // hit water surface

        weight *= attenuation_rate_water; 
        TraceShadowRay(hit_point, dir, inv_dir, weight);
        return;
    }

}


void SPH::RayCasting(double *org, double *dir, double *inv_dir, int depth, double *color)
{
    if(depth > 1){
        color[0] = color[1] = color[2] = 1.0;
        return;
    }

    double hit_point[3], normal[3];
    
    int nearest_object_index = FindNearestObject(org, dir, inv_dir, hit_point, normal);


    if(nearest_object_index == -1){
        color[0] = color[1] = color[2] = 1.0;
        return;
    }else if(nearest_object_index >= 0 && nearest_object_index <= 3){

        ComputeDirectLighting(hit_point, normal, dir, nearest_object_index, depth, color);
        return;
    
    }else if(nearest_object_index == 4){

        ComputeDirectLighting(hit_point, normal, dir, nearest_object_index, depth, color);      

        double back_color[3];
        RayCasting(hit_point, dir, inv_dir, depth, back_color);

        for(int i = 0; i < 3; i++) color[i] = color[i]*0.2 + back_color[i]*0.8;
        
        return;
    
    }else if(nearest_object_index == 5){

        ComputeDirectLighting(hit_point, normal, dir, nearest_object_index, depth, color);      
        return;
    
    }else if(nearest_object_index == 6){

        ComputeDirectLighting(hit_point, normal, dir, nearest_object_index, depth, color);      
        return;
    
    }else if(nearest_object_index == 100){

        ComputeDirectLighting(hit_point, normal, dir, nearest_object_index, depth, color);

        double refraction_vector[3];
        ComputeRefractionVector(normal, dir, refraction_vector);

        for(int i = 0; i < 3; i++) dir[i] = refraction_vector[i];
        for(int i = 0; i < 3; i++) inv_dir[i] = 1.0/dir[i];

        double back_color[3];
        RayCasting(hit_point, dir, inv_dir, depth, back_color);

        for(int i = 0; i < 3; i++) color[i] = color[i]*0.1 + back_color[i]*0.9;

        return;
    }

}


void SPH::RenderUsingRaycasting(double *eye, int window_width, int window_height)
{
    static float *image;
    bool init = true;

    if(init){
        image = new float [window_width*window_height*3];
        init = false;
    }

    double rotated_org[3] = { 0.0, 0.0, 0.0, };

    float rot_matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, rot_matrix);

    for(int i = 0; i < 3; i++){
        rotated_org[0] += rot_matrix[i]   * eye[i];
        rotated_org[1] += rot_matrix[4+i] * eye[i];
        rotated_org[2] += rot_matrix[8+i] * eye[i];
    }


    double target_origin[2] = { -sin(20.0/180.0*PI), -sin(20.0/180.0*PI)  };

    static double deltaW = (2.0*sin(20.0/180.0*PI)) / window_width;
    static double deltaH = (2.0*sin(20.0/180.0*PI)) / window_height;


    for (int i = 0; i < window_height; i++){

        for (int j = 0; j < window_width; j++) { 

            double color[3] = { 0.0, 0.0, 0.0, };

            for(int m = 0; m < 1; m++){
                for(int n = 0; n < 1; n++){

                    double dir[3] = { target_origin[0]+deltaW*(j-0.25*m), target_origin[1]+deltaH*(i-0.25*n), -1.0 };

                    double rotated_dir[3] = { 0.0, 0.0, 0.0, };
                    for(int k = 0; k < 3; k++){
                        rotated_dir[0] += rot_matrix[k]   * dir[k];
                        rotated_dir[1] += rot_matrix[4+k] * dir[k];
                        rotated_dir[2] += rot_matrix[8+k] * dir[k];
                    }
                    Normalize(rotated_dir);

                    double inv_rotated_dir[3];
                    for(int k = 0; k < 3; k++) inv_rotated_dir[k] = 1.0 / rotated_dir[k];


                    RayCasting(rotated_org, rotated_dir, inv_rotated_dir, 0, color);  

                    image[(i*window_width+j)*3]    += color[0];// * 0.0625;
                    image[(i*window_width+j)*3 +1] += color[1];// * 0.0625;
                    image[(i*window_width+j)*3 +2] += color[2];// * 0.0625;
                }
            }

        }

    }


    glDrawPixels(window_width, window_height, GL_RGB, GL_FLOAT, image);

    wchar_t filename[100];
    static int file_id_base = 0;

    swprintf(filename, 100, L"img/sample_%d.png", file_id_base++);

    SaveImage(filename, window_width, window_height, image);
}
