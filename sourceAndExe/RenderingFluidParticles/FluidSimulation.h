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

#include "octree.h"

struct WaterParticle {
    double position[3], velocity[3], acceleration[3];
    double inv_density;
    double pressure;
    double gravity_force;

    float color[3];

    WaterParticle(double *p_in, double *v_in){
        for(int i = 0; i < 3; i++) position[i]     = p_in[i];
        for(int i = 0; i < 3; i++) velocity[i]     = v_in[i];
        for(int i = 0; i < 3; i++) acceleration[i] = 0.0;
    }
};

struct Grid {
    double bbmin[3], bbmax[3];
    bool   isObstacle;

    vector<WaterParticle> water_particles;

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

    double surface_distance;

    CellVertex(){};
};



class SPH {

    Octree octree;

    int    FindNearestObject(double *org, double *dir, double *inv_dir, double *hit_point, double *normal);
    void   TraceShadowRay(double *org, double *dir, double *inv_dir, double &weight);
    void   ComputeDirectLighting(double *point, double *normal, double *direction, int object_index, double depth, double *color);
    void   RayCasting(double *org, double *dir, double *inv_dir, int depth, double *color);

public:
    SPH(){ }
    void InitializeGrids();
    void ProcessStep();
    void RenderParticles();
    void OutputParticleLocation();
    void RenderUsingRaycasting(double *eye, int window_width, int window_height);
};

typedef vector<WaterParticle>::iterator  WaterParticleIter;



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
extern void ComputeRefractionVector(double *normal, double *vector_in, double *reraction_vector);
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

#define EPSILON 1.0e-6



