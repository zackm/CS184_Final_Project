#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>
#include <stack>
#include <utility>
#include <GL/glut.h>
using namespace std;

#include <windows.h>
#include <gdiplus.h>
using namespace Gdiplus;
#pragma comment (lib,"Gdiplus.lib")



struct Particle {
    double position[3], velocity[3], acceleration[3];
    double inv_density;
    double pressure;
    double gravity_force;
    double buoyant_force;

    bool  is_outside;


    Particle(double *p_in, double *v_in){
        for(int i = 0; i < 3; i++) position[i]     = p_in[i];
        for(int i = 0; i < 3; i++) velocity[i]     = v_in[i];
        for(int i = 0; i < 3; i++) acceleration[i] = 0.0;

        is_outside = false;
    }
};

struct Grid {
    double bbmin[3], bbmax[3];

    vector<Particle> water_particles;
    vector<Particle> air_particles;

    Grid(){};
    Grid(double *bbmin_in, double *bbmax_in){
        for(int i = 0; i < 3; i++) bbmin[i] = bbmin_in[i];
        for(int i = 0; i < 3; i++) bbmax[i] = bbmax_in[i];
    }
};

struct CellVertex {
    double color_field;
    double color_field_normal[3];
    double color_field_normal_norm;
    double color_field_laplacian;

    CellVertex(){};
};

//#define GRID_Y 1000

class SPH {
    Grid       grid[28][20][5];
    CellVertex cell_vertex[29][21][6];
    int  n_particles;

    void UpdateDensity(int particle_id, int id_x, int id_y, int id_z, bool is_water_particle);
    void UpdatePartcle(int particle_id, int id_x, int id_y, int id_z, bool is_water_particle, vector<Particle> &new_air_particles);
    void ProcessParticleSolidCollision(double *new_position, double *pos, double *vel);
    void UpdateGrid(int id_x, int id_y, int id_z);

    void ComputeDensity(int particle_id, int id_x, int id_y, int id_z);
    void ComputeColorField(CellVertex &cell_vertex, int id_x, int id_y, int id_z);

public:
    SPH(){ n_particles = 0; }
    void InitializeGrids();
    void ProcessStep();
    void RenderParticles();
    void OutputParticleLocation();
};

typedef vector<Particle>::iterator  WaterParticleIter;



extern void GLInit();
extern void CrossProduct( double *a, double *b, double *c );
extern void Normalize(double *a);
extern double DotProduct(double *a, double *b);
extern double DotProduct4D(double *a, double *b);
extern void GetEdgeVector(double *v1, double *v2, double *edge_vector);
extern void GetArea(double *normal, double &area);
extern double GetDistance(double *a, double *b);
extern double GetLength(double *a);
extern void Swap(double &a, double &b);
extern double Max(double a, double b);
extern double Min(double a, double b);
extern bool SolveLinearSystem(double (*matrix)[4], double *rhs, double *solution);
extern void RotateAroundAxis(double *rotation_axis, double theta, double *vec_in, double *vec_out);
extern void ComputeReflectionVector(double *normal, double *vector_in, double *reflection_vector);
extern bool ComputeRayBoxIntersection(double *org, double *dir, double *inv_dir, double *bbmin, double *bbmax, 
                                      double &nearest_t, double *intersection, double *normal);
extern double random0to1();
extern float* LoadPngData(WCHAR *filename, unsigned int &img_width, unsigned int &img_height);
extern void SaveImage(wchar_t *filename, int w, int h, float *buffer);


extern double Wpoly6(double r_square);
extern void   WspikyGrad(double *r, double r_square, double *gradient);
extern double WviscosityLaplacian(double r_square);
extern void   Wpoly6Grad(double *r, double r_square, double *gradient);
extern double Wpoly6Laplacian(double r_square);

extern GLint window_width, window_height;
extern float *png_data;
extern unsigned int   png_width, png_height;


#define PI    3.14159265
#define TwoPI 6.28318531
#define InvPI 0.318309886

#define H 0.045
#define Rho0_w 1000 // (kg/m^3) water particle rest density 
#define Rho0_a 10
#define WaterParticleMass 0.012 // (kg)
#define AirParticleMass   0.0012 // (kg)
#define Sigma_s 0.6  // surface tension coefficient

#define DeltaT 0.001 // time step

extern double K;
extern double Mu_w, Mu_a;
extern bool is_air_trap;

#define EPSILON 1.0e-6



