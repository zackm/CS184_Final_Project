
#include "FluidSimulation.h"

double K  = 20;
double Mu_w = 5;
double Mu_a = 1;

bool is_air_trap = true;


void SPH::UpdateDensity(int particle_id, int id_x, int id_y, int id_z, bool is_water_particle)
{
    double *pos;

    if(is_water_particle) pos = grid[id_x][id_y][id_z].water_particles[particle_id].position;
    else                  pos = grid[id_x][id_y][id_z].air_particles[particle_id].position;

    double sumW_water = 0.0, sumW_air = 0.0;

    for(int x = -1; x <= 1; x++){
        if(id_x+x < 0)   continue;
        if(id_x+x >= 28) break;

        for(int y = -1; y <= 1; y++){
            if(id_y+y < 0)   continue;
            if(id_y+y >= 20) break;

            for(int z = -1; z <= 1; z++){
                if(id_z+z < 0)  continue;
                if(id_z+z >= 4) break;

                for(unsigned int i = 0; i < grid[id_x+x][id_y+y][id_z+z].water_particles.size(); i++){

                    double *neighbor_pos = grid[id_x+x][id_y+y][id_z+z].water_particles[i].position;
                    double r[3] = { pos[0]-neighbor_pos[0], pos[1]-neighbor_pos[1], pos[2]-neighbor_pos[2] };

                    double r_square = DotProduct(r, r);

                    if(r_square <= H*H) sumW_water += Wpoly6(r_square);
                }

                for(unsigned int i = 0; i < grid[id_x+x][id_y+y][id_z+z].air_particles.size(); i++){

                    double *neighbor_pos = grid[id_x+x][id_y+y][id_z+z].air_particles[i].position;
                    double r[3] = { pos[0]-neighbor_pos[0], pos[1]-neighbor_pos[1], pos[2]-neighbor_pos[2] };

                    double r_square = DotProduct(r, r);

                    if(r_square <= H*H) sumW_air += Wpoly6(r_square);
                }

            } // for(int z = -1; z <= 1; z++
        } // for(int y = -1; y <= 1; y++){
    } // for(int x = -1; x <= 1; x++){


    double density = WaterParticleMass*sumW_water + AirParticleMass*sumW_air;  // "sumW" corresnponds to volume??

    if(is_water_particle){
        grid[id_x][id_y][id_z].water_particles[particle_id].inv_density   = 1.0 / density; 
        grid[id_x][id_y][id_z].water_particles[particle_id].pressure      = K*(density - Rho0_w); 
        grid[id_x][id_y][id_z].water_particles[particle_id].gravity_force = -density*9.80665;
        grid[id_x][id_y][id_z].water_particles[particle_id].buoyant_force = 0.0;
    }else{
        grid[id_x][id_y][id_z].air_particles[particle_id].inv_density   = 1.0 / density; 
        grid[id_x][id_y][id_z].air_particles[particle_id].pressure      = K*(density - Rho0_a); 
        grid[id_x][id_y][id_z].air_particles[particle_id].gravity_force = -density*9.80665;
        grid[id_x][id_y][id_z].air_particles[particle_id].buoyant_force = 13.0*(density-Rho0_a)*9.80665;
    }

}


void SPH::UpdatePartcle(int particle_id, int id_x, int id_y, int id_z, bool is_water_particle, vector<Particle> &new_air_particles)
{
    Particle *pt;
    double viscosity;


    if(is_water_particle){
        pt = &(grid[id_x][id_y][id_z].water_particles[particle_id]);
        viscosity = Mu_w;
    }else{
        pt = &(grid[id_x][id_y][id_z].air_particles[particle_id]);
        viscosity = Mu_a;
    }

    double f_pressure_w[3] = { 0.0, 0.0, 0.0 }, f_viscosity_w[3] = { 0.0, 0.0, 0.0 };
    double f_pressure_a[3] = { 0.0, 0.0, 0.0 }, f_viscosity_a[3] = { 0.0, 0.0, 0.0 };

    double color_field_normal_w[3] = { 0.0, 0.0, 0.0 }, color_field_laplacian_w = 0.0;
    double color_field_normal_a[3] = { 0.0, 0.0, 0.0 }, color_field_laplacian_a = 0.0;


    for(int x = -1; x <= 1; x++){
        if(id_x+x < 0)   continue;
        if(id_x+x >= 28) break;

        for(int y = -1; y <= 1; y++){
            if(id_y+y < 0)   continue;
            if(id_y+y >= 20) break;

            for(int z = -1; z <= 1; z++){
                if(id_z+z < 0)  continue;
                if(id_z+z >= 4) break;

                for(unsigned int i = 0; i < grid[id_x+x][id_y+y][id_z+z].water_particles.size(); i++){

                    Particle *neightbor_pt = &(grid[id_x+x][id_y+y][id_z+z].water_particles[i]);

                    double r[3] = { pt->position[0] - neightbor_pt->position[0], 
                                    pt->position[1] - neightbor_pt->position[1], 
                                    pt->position[2] - neightbor_pt->position[2] };
                    double r_square = DotProduct(r, r);

                    if(r_square <= H*H){

                        if(r_square > 0.0){
                            double gradient[3];
                            Wpoly6Grad(r, r_square, gradient);

                            double weight_p = (pt->pressure + neightbor_pt->pressure)*0.5 * neightbor_pt->inv_density;

                            for(int j = 0; j < 3; j++) f_pressure_w[j] += weight_p * gradient[j];

                            for(int j = 0; j < 3; j++) color_field_normal_w[j] += neightbor_pt->inv_density * gradient[j];
                        }

                        double weight_v = neightbor_pt->inv_density * WviscosityLaplacian(r_square);

                        for(int j = 0; j < 3; j++) f_viscosity_w[j] += (neightbor_pt->velocity[j] - pt->velocity[j]) * weight_v;                     

                        color_field_laplacian_w += neightbor_pt->inv_density * Wpoly6Laplacian(r_square);
                    } // if(r_square <= H*H){

                }

                for(unsigned int i = 0; i < grid[id_x+x][id_y+y][id_z+z].air_particles.size(); i++){

                    Particle *neightbor_pt = &(grid[id_x+x][id_y+y][id_z+z].air_particles[i]);

                    double r[3] = { pt->position[0] - neightbor_pt->position[0], 
                                    pt->position[1] - neightbor_pt->position[1], 
                                    pt->position[2] - neightbor_pt->position[2] };
                    double r_square = DotProduct(r, r);

                    if(r_square <= H*H){

                        if(r_square > 0.0){
                            double gradient[3];
                            Wpoly6Grad(r, r_square, gradient);

                            double weight_p = (pt->pressure + neightbor_pt->pressure)*0.5 * neightbor_pt->inv_density;

                            for(int j = 0; j < 3; j++) f_pressure_a[j] += weight_p * gradient[j];
                            
                            for(int j = 0; j < 3; j++) color_field_normal_a[j] += neightbor_pt->inv_density * gradient[j];
                        }

                        double weight_v = neightbor_pt->inv_density * WviscosityLaplacian(r_square);

                        for(int j = 0; j < 3; j++) f_viscosity_a[j] += (neightbor_pt->velocity[j]-pt->velocity[j]) * weight_v;

                        color_field_laplacian_a += neightbor_pt->inv_density * Wpoly6Laplacian(r_square);
                    } // if(r_square <= H*H){

                }

            } // for(int z = -1; z <= 1; z++
        } // for(int y = -1; y <= 1; y++){
    } // for(int x = -1; x <= 1; x++){


    double f_pressure[3], f_viscosity[3];
    double f_surface[3] = { 0.0, 0.0, 0.0 };
    double color_field_normal[3], color_field_laplacian;

    for(int j = 0; j < 3; j++) f_pressure[j]  = -WaterParticleMass*f_pressure_w[j]  -AirParticleMass*f_pressure_a[j];
    for(int j = 0; j < 3; j++) f_viscosity[j] = (viscosity+Mu_w)*0.5*WaterParticleMass*f_viscosity_w[j] + (viscosity+Mu_a)*0.5*AirParticleMass*f_viscosity_a[j];
    for(int j = 0; j < 3; j++) color_field_normal[j] = color_field_normal_w[j]*WaterParticleMass + color_field_normal_a[j]*AirParticleMass;
    color_field_laplacian = color_field_laplacian_w*WaterParticleMass + color_field_laplacian_a*AirParticleMass;

    double normal_norm = sqrt( DotProduct(color_field_normal, color_field_normal) );


    pt->acceleration[0] = (f_pressure[0]+f_viscosity[0]+f_surface[0]) * pt->inv_density;
    pt->acceleration[1] = (f_pressure[1]+f_viscosity[1]+f_surface[1]+pt->gravity_force+pt->buoyant_force) * pt->inv_density;
    pt->acceleration[2] = (f_pressure[2]+f_viscosity[2]+f_surface[2]) * pt->inv_density;

    // update position 
    double prev_pos[3];

    for(int j = 0; j < 3; j++) prev_pos[j] = pt->position[j];
        
    for(int j = 0; j < 3; j++) pt->position[j] += (pt->velocity[j] + pt->acceleration[j]*0.5*DeltaT) * DeltaT;
    for(int j = 0; j < 3; j++) pt->velocity[j] += pt->acceleration[j] * DeltaT;



    if(pt->position[1] > 0.899 && (pt->position[0] <= -0.45 || pt->position[0] >= 0.45 || !is_air_trap)){
        if(is_water_particle)
            grid[id_x][id_y][id_z].water_particles[particle_id].is_outside = true;
        else
            grid[id_x][id_y][id_z].air_particles[particle_id].is_outside = true;
    }


    // process collision between particle and walls and confirm final new position and velocity
    ProcessParticleSolidCollision(prev_pos, pt->position, pt->velocity);

    if(is_air_trap && is_water_particle && normal_norm > 21.0){ 
        // generate air particle
        Normalize(color_field_normal);

        double d = 0.05;
        double init_position[3] = { pt->position[0]-color_field_normal[0]*d, 
                                    pt->position[1]-color_field_normal[1]*d, 
                                    pt->position[2]-color_field_normal[2]*d };
        double init_velocity[3] = { pt->velocity[0], pt->velocity[1], pt->velocity[2] };

        ProcessParticleSolidCollision(pt->position, init_position, init_velocity);

        if(-color_field_normal[1] > 0.99 && id_x >=8 && id_x < 19 && id_y >=4 && id_y < 16){
            new_air_particles.push_back( Particle(init_position, init_velocity) );
        }

    }

    if( !is_water_particle && (id_x < 8 || id_x >= 20 || id_y < 4 || id_y > 16) )
         grid[id_x][id_y][id_z].air_particles[particle_id].is_outside = true;

}

void SPH::ProcessParticleSolidCollision(double *prev_pos, double *pos, double *vel)
{
    double dir[3], inv_dir[3];

    for(int i = 0; i < 3; i++) dir[i]     = pos[i]-prev_pos[i];
    for(int i = 0; i < 3; i++) inv_dir[i] = 1.0 / dir[i];

    double nearest_t = 1.0e6;
    int    xyz = -1;

    double bbmin[5][3] = { {-0.63,  0.0, 0.0}, {-0.45, 0.18, 0.0, },  {-0.09, 0.0, 0.0, }, { 0.27, 0.18, 0.0, }, {-0.27, 0.72, 0.0, } };
    double bbmax[5][3] = { { 0.63, 0.90, 0.18}, {-0.27, 0.90,  0.18, }, { 0.09, 0.36,  0.18, }, { 0.45, 0.90,  0.18, }, { 0.27, 0.90,  0.18, },} ;


    for(int i = 0; i < 5; i++){
        if(!is_air_trap && i == 4) break;

        double t, intersection[3], normal[3];

        if( ComputeRayBoxIntersection(prev_pos, dir, inv_dir, bbmin[i], bbmax[i], t, intersection, normal) ){

            if(t > 0.0 && t < 1.0+EPSILON && t < nearest_t){ 
                nearest_t = t;

                for(int j = 0; j < 3; j++) pos[j] = prev_pos[j] + dir[j]*(t-EPSILON);
                for(int j = 0; j < 3; j++){
                    if(fabs(normal[j]) > 1.0e-6){
                        xyz = j; 
                    }
                }
            }
        }
    }

    if(nearest_t > 0.0 && nearest_t < 1.0+EPSILON){ // collision occur
        vel[xyz] = -vel[xyz]*0.5;
    }

}




void SPH::UpdateGrid(int id_x, int id_y, int id_z)
{
    static double inv_x_step = 1.0/0.045;
    static double inv_y_step = 1.0/0.045;
    static double inv_z_step = 1.0/0.045;


    for(int i = 0; i < (int)grid[id_x][id_y][id_z].water_particles.size(); i++){

        if(grid[id_x][id_y][id_z].water_particles[i].is_outside){
            // delete particle
            grid[id_x][id_y][id_z].water_particles[i] = grid[id_x][id_y][id_z].water_particles.back();
            grid[id_x][id_y][id_z].water_particles.pop_back();
            i--;
            continue;
        }

        double *pos = grid[id_x][id_y][id_z].water_particles[i].position;

        int current_id_x = (int)floor( (pos[0]+0.63)*inv_x_step );
        int current_id_y = (int)floor( (pos[1])     *inv_y_step );
        int current_id_z = (int)floor( (pos[2])     *inv_z_step );


        if(current_id_x < 0)        current_id_x = 0;
        else if(current_id_x >= 28) current_id_x = 27;
        if(current_id_y < 0)        current_id_y = 0;
        else if(current_id_y >= 20) current_id_y = 19;
        if(current_id_z < 0)        current_id_z = 0;
        else if(current_id_z >= 4 ) current_id_z = 3;


        if(id_x != current_id_x || id_y != current_id_y || id_z != current_id_z){
            // move particle to new grid
            grid[current_id_x][current_id_y][current_id_z].water_particles.push_back( grid[id_x][id_y][id_z].water_particles[i] );

            // delete the particle from this grid
            grid[id_x][id_y][id_z].water_particles[i] = grid[id_x][id_y][id_z].water_particles.back();
            grid[id_x][id_y][id_z].water_particles.pop_back();
            i--;
        }

    }


    for(int i = 0; i < (int)grid[id_x][id_y][id_z].air_particles.size(); i++){

        if(grid[id_x][id_y][id_z].air_particles[i].is_outside){
            // delete particle
            grid[id_x][id_y][id_z].air_particles[i] = grid[id_x][id_y][id_z].air_particles.back();
            grid[id_x][id_y][id_z].air_particles.pop_back();
            i--;
            continue;
        }

        double *pos = grid[id_x][id_y][id_z].air_particles[i].position;

        int current_id_x = (int)floor( (pos[0]+0.63)*inv_x_step );
        int current_id_y = (int)floor( (pos[1])     *inv_y_step );
        int current_id_z = (int)floor( (pos[2])     *inv_z_step );


        if(current_id_x < 0)        current_id_x = 0;
        else if(current_id_x >= 28) current_id_x = 27;
        if(current_id_y < 0)        current_id_y = 0;
        else if(current_id_y >= 20) current_id_y = 19;
        if(current_id_z < 0)        current_id_z = 0;
        else if(current_id_z >= 4 ) current_id_z = 3;

        if(id_x != current_id_x || id_y != current_id_y || id_z != current_id_z){
            // move particle to new grid
            grid[current_id_x][current_id_y][current_id_z].air_particles.push_back( grid[id_x][id_y][id_z].air_particles[i] );

            grid[id_x][id_y][id_z].air_particles[i] = grid[id_x][id_y][id_z].air_particles.back();
            grid[id_x][id_y][id_z].air_particles.pop_back();
            i--;
        }

    }

}


void SPH::ProcessStep()
{

    if(n_particles < 22000){
        for(int i = 0; i < 2; i++){
            double init_position[3] = { -0.36, 0.025, 0.09 };
            double init_velocity[3] = { random0to1(), random0to1(), random0to1() };

            grid[6][0][3].water_particles.push_back( Particle(init_position, init_velocity) );
            n_particles++;
        }
    }


    for(int x = 0; x < 28; x++){
        for(int y = 0; y < 20; y++){
            for(int z = 0; z < 4; z++){

                for(unsigned int i = 0; i < grid[x][y][z].water_particles.size(); i++){
                    UpdateDensity(i, x, y, z, true);
                }
                for(unsigned int i = 0; i < grid[x][y][z].air_particles.size(); i++){
                    UpdateDensity(i, x, y, z, false);
                }

            }
        }
    }

    vector<Particle> new_air_particles[28][20][4];
    

    for(int x = 0; x < 28; x++){
        for(int y = 0; y < 20; y++){
            for(int z = 0; z < 4; z++){

                for(unsigned int i = 0; i < grid[x][y][z].water_particles.size(); i++){
                    UpdatePartcle(i, x, y, z, true, new_air_particles[x][y][z]);
                }
                for(unsigned int i = 0; i < grid[x][y][z].air_particles.size(); i++){
                    UpdatePartcle(i, x, y, z, false, new_air_particles[x][y][z]);
                }

            }
        }
    }


    for(int x = 0; x < 28; x++){
        for(int y = 0; y < 20; y++){
            for(int z = 0; z < 4; z++){
                grid[x][y][z].air_particles.insert(grid[x][y][z].air_particles.end(), 
                                                   new_air_particles[x][y][z].begin(), new_air_particles[x][y][z].end());
            }
        }
    }


    for(int x = 0; x < 28; x++){
        for(int y = 0; y < 20; y++){
            for(int z = 0; z < 4; z++){
                UpdateGrid(x, y, z);
            }
        }
    }

}

void SPH::RenderParticles()
{
    glPointSize(4);
    glEnable(GL_LIGHTING);

    for(int x = 0; x < 28; x++){
        for(int y = 0; y < 20; y++){
            for(int z = 0; z < 4; z++){

                float water_particle_color[3] = {0,0,1}, air_particle_color[3] = {1,0.6,0};

                for(unsigned int i = 0; i < grid[x][y][z].water_particles.size(); i++){
                    double *pos = grid[x][y][z].water_particles[i].position;

                    glMaterialfv(GL_FRONT, GL_DIFFUSE, water_particle_color); 

                    glPushMatrix();
                    glTranslated(pos[0], pos[1], pos[2]);
                    glutSolidSphere(0.006, 4, 4);
                    glPopMatrix();
                }

                for(unsigned int i = 0; i < grid[x][y][z].air_particles.size(); i++){

                    double *pos = grid[x][y][z].air_particles[i].position;

                    glMaterialfv(GL_FRONT, GL_DIFFUSE, air_particle_color); 

                    glPushMatrix();
                    glTranslated(pos[0], pos[1], pos[2]);
                    glutSolidSphere(0.006, 4, 4);
                    glPopMatrix();
                }

            }
        }
    }

}

void SPH::OutputParticleLocation()
{
    char filename_out[100];

    static int file_base = 0;

    // stop simulation
    if(file_base == 400) exit(0);

    sprintf_s(filename_out, 100, "data/pt_data_%d", file_base++);
    
    cerr << "saving: " << file_base << " / 400" << endl;

    ofstream particle_location_data(filename_out, ios::binary);

    particle_location_data.write((char*)&n_particles, sizeof(int));

    for(int x = 0; x < 28; x++){
        for(int y = 0; y < 20; y++){
            for(int z = 0; z < 4; z++){

                for(unsigned int i = 0; i < grid[x][y][z].water_particles.size(); i++){

                    double *pos = grid[x][y][z].water_particles[i].position;

                    particle_location_data.write((char*)&pos[0], sizeof(double));
                    particle_location_data.write((char*)&pos[1], sizeof(double));
                    particle_location_data.write((char*)&pos[2], sizeof(double));
                }

            }
        }
    }

    particle_location_data.close();

}


void SPH::InitializeGrids()
{

    for(int x = 0; x < 28; x++){
        for(int y = 0; y < 20; y++){
            for(int z = 0; z < 4; z++){
                double bbmin[3] = {     x*0.045-0.63,     y*0.045,     z*0.045 };
                double bbmax[3] = { (x+1)*0.045-0.63, (y+1)*0.045, (z+1)*0.045 };

                for(int i = 0; i < 3; i++) grid[x][y][z].bbmin[i] = bbmin[i];
                for(int i = 0; i < 3; i++) grid[x][y][z].bbmax[i] = bbmax[i]; 
            }
        }
    }

}
