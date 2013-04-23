
#include "FluidSimulation.h"



void SPH::InitializeGrids()
{
    vector<WaterParticle> water_particles;

    char filename_in[100];

    static int file_base = 0;
    sprintf_s(filename_in, 100, "data/pt_data_%d", file_base++);

    cerr << filename_in << endl;

    ifstream fin(filename_in, ios::binary);

    if(!fin.is_open()) exit(0);

    int n_particles = 0;
    fin.read((char*)&n_particles, sizeof(int));

    for(int i = 0; i < n_particles; i++){
        double pos[3], dummy[3];
        fin.read((char*)&pos, sizeof(double)*3);

        if(fin.eof()) break;

        water_particles.push_back( WaterParticle(pos, dummy) );
    }


    octree.DeleteOctree();
    octree.InitOctree(water_particles);

    cerr << "n_particles " << n_particles << endl;
}
