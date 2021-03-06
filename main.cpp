#include "main.h"

#include "FreeImage.h"

#include <sstream>
#include <fstream>
#include <cstring>

using namespace std;

/*******************
* GLOBAL VARIABLES *
*******************/
Container CONTAINER(Vec3(2,2,2),Vec3(0,0,0));
vector<Particle*> PARTICLES;//particles that we do SPH on.
vector<Triangle*> TRIANGLES;//triangles made from marching cubes to render

const float TIMESTEP = .01;//time elapsed between iterations
const float LIFETIME = 100.0f;
float CURRENT_TIME = 0.0f;
int NUM_PARTICLES = 0;
Vec3 GRAVITY(0,-9.8f,0);
const float MASS = .02f;//could set it to any number really.
const float IDEAL_DENSITY = 1000.0f;
const float STIFFNESS = 3.0f;//for pressure difference
float VISCOSITY = 3.5f;
float SURFACE_TENSION = .07f;
const float TENSION_THRESHOLD = 7.0f;

float CUBE_TOL = .01f;//either grid size or tolerance for adaptive cubes, reciprocal must be an integer for now. .025f ok, .01 for high quality
const float DENSITY_TOL = 100.0f;//also used for marching grid, for density of the particles
const int WIDTH = 20;//floor((CONTAINER.max.x-CONTAINER.min.x)/CUBE_TOL);
const int HEIGHT = 20;//floor((CONTAINER.max.y-CONTAINER.min.y)/CUBE_TOL);
const int DEPTH = 20;//floor((CONTAINER.max.z-CONTAINER.min.z)/CUBE_TOL);


Neighbor NEIGHBOR; //neighbor object used for calculations
const float H = .0625; // .0625 works well, .05 good too
const float SUPPORT_RADIUS = .125; // .125 works well, .1 good too

bool RENDERING_TRIANGLES = false;
bool RENDERING_BLOB = false;
bool RENDERING_WIREFRAME = false;

const float PI = 3.1415926;
const float DRAW_RADIUS = .01f;

bool OUTPUT_IMAGE = false; // boolean to contstantly output opengl frames for movie
bool OUTPUT_SINGLE_IMAGE = false; // boolean to output a single opengl frame, continue execution
bool RAYTRACE_MOVIE = false; // boolean to constantly output raytracing frames for movie
int IMAGE_COUNTER = 0; // keep track of frame number for outputting images and obj files
int IMAGE_DELAY = 0; // number of frames to wait before outputting images
int PIC_WIDTH = 600; // display window and output image dimensions
int PIC_HEIGHT = 600;
bool HIDE = false;

bool SHOW_FPS = false; // print out FPS to console
int SETUP_SCENE = 0; // choice of preset scene, from command line input
bool PAUSE = false; // run time, or pause

Vec3 normal_at_point(Vec3);

//marching cubes table data
int edgeTable[256]={
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };
char triTable[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

/* 
Output triangle mesh to OBJ file.
*/
void output_obj() {
    // place particles in box_neighbor
    NEIGHBOR.place_particles(PARTICLES, SUPPORT_RADIUS, CONTAINER, NUM_PARTICLES, true);
    
	// open file
	std::stringstream ss1;
	ss1 << IMAGE_COUNTER;
	std::string s1(ss1.str());
	string save_name = std::string("Multi_Trace/input_obj/fluid")+s1+".obj";

	ofstream output_file;
	output_file.open (save_name.c_str());
	output_file << "# OBJ File created by Tyler Brabham and Zack Mayeda\n";
	output_file << "# UC Berkeley CS184 Spring 2013 Final Project\n\n";

	if(RENDERING_TRIANGLES){
		vector<Vec3> v,vn;

		for (int i = 0; i < TRIANGLES.size(); i++) {
			Triangle *t = TRIANGLES[i];
			Vec3 v1(t->a.x,t->a.y,t->a.z);
			Vec3 v2(t->b.x,t->b.y,t->b.z);
			Vec3 v3(t->c.x,t->c.y,t->c.z);
			Vec3 vn_1(t->a_normal.x,t->a_normal.y,t->a_normal.z);
			Vec3 vn_2(t->b_normal.x,t->b_normal.y,t->b_normal.z);
			Vec3 vn_3(t->c_normal.x,t->c_normal.y,t->c_normal.z);
			v.push_back(v1); v.push_back(v2); v.push_back(v3);
			vn.push_back(vn_1); vn.push_back(vn_2); vn.push_back(vn_3);
		}
		// write all vertices, v x y z
		for (int i = 0; i < v.size(); i++) {
			Vec3 v_temp = v[i];
			output_file<<"v "<< v_temp.x<<" "<<v_temp.y<< " "<<v_temp.z<<endl;
		}
		// write all vertex normals, vn x y z
		for (int i = 0; i < vn.size(); i++) {
			Vec3 vn_temp = vn[i];
			output_file<<"vn "<< vn_temp.x<<" "<<vn_temp.y<< " "<<vn_temp.z<<endl;
		}
	}else{
		//we are rendering blobs, so we will make an obj point file.
		Particle* temp_part;
		for (int i = 0; i < NUM_PARTICLES; i++) {
			temp_part = PARTICLES[i];
			output_file<<"part "<< temp_part->position.x<<" "<<temp_part->position.y<<" "<<temp_part->position.z
				<<" "<<temp_part->mass<<' '<<temp_part->density<<endl;
		}

	}

	output_file.close();
    
	Raytracer r(NEIGHBOR);

    std::stringstream ss;
    ss << IMAGE_COUNTER;
    std::string s(ss.str());
    string out_name = std::string("Multi_Trace/output_pics/fluid")+s+".png";
	r.ray_trace_start(save_name,out_name,1200,1200);
    IMAGE_COUNTER++;
}

/* 
Keyboard interactions.
*/
void keyPressed(unsigned char key, int x, int y) {
	switch(key) {
	case ' ':
        // Pause program
            PAUSE = !PAUSE;
		break;
    case 'q':
        // Exit Program
        exit(0);
        break;
	case 'r':
        // Output a single obj file, raytrace it, and quit
		output_obj();
		exit(0);
		break;
	case 'p':
        // Output a single opengl frame and continue
		OUTPUT_SINGLE_IMAGE = true;
		break;
	case 'i':
        // Toggle isosurface rendering
		RENDERING_TRIANGLES = !RENDERING_TRIANGLES;
		break;
	case 'f':
        // Print out FPS to console once per keypress
		SHOW_FPS = true;
		break;
	}
}

/*
Linearly interpolate the position where an isosurface cuts
an edge between two vertices, each with their own scalar value.
This code comes from Matt Fisher at http://graphics.stanford.edu/~mdfisher/MarchingCubes.html
It is cited in our final project write up.
*/
Vec3 VertexInterp(Vec3 *p1,Vec3* p2,float valp1,float valp2){
	return ((*p1) + ((*p2) - (*p1))*(-valp1 / (valp2 - valp1)));
}

/*
This code generates the list of triangles at the gridcell. It comes from Matt Fisher, and
is again cited in the final project write up.
*/
vector<Triangle*> polygonise(GRIDCELL &Grid, int &NewVertexCount, vector<Vec3> vertices){
	int TriangleCount;
	int CubeIndex;
	vector<Triangle*> new_triangles;
	vector<Vec3> VertexList;
	vector<Vec3> NewVertexList;
	char LocalRemap[12];

	VertexList.resize(12);
	NewVertexList.resize(12);

	//Determine the index into the edge table which
	//tells us which vertices are inside of the surface
	CubeIndex = 0;
	if (Grid.val[0] < 0.0f) CubeIndex |= 1;
	if (Grid.val[1] < 0.0f) CubeIndex |= 2;
	if (Grid.val[2] < 0.0f) CubeIndex |= 4;
	if (Grid.val[3] < 0.0f) CubeIndex |= 8;
	if (Grid.val[4] < 0.0f) CubeIndex |= 16;
	if (Grid.val[5] < 0.0f) CubeIndex |= 32;
	if (Grid.val[6] < 0.0f) CubeIndex |= 64;
	if (Grid.val[7] < 0.0f) CubeIndex |= 128;

	//Cube is entirely in/out of the surface
	if (edgeTable[CubeIndex] == 0)
		return new_triangles;

	//Find the vertices where the surface intersects the cube
	if (edgeTable[CubeIndex] & 1)
		VertexList[0] =
		VertexInterp(Grid.p[0],Grid.p[1],Grid.val[0],Grid.val[1]);
	if (edgeTable[CubeIndex] & 2)
		VertexList[1] =
		VertexInterp(Grid.p[1],Grid.p[2],Grid.val[1],Grid.val[2]);
	if (edgeTable[CubeIndex] & 4)
		VertexList[2] =
		VertexInterp(Grid.p[2],Grid.p[3],Grid.val[2],Grid.val[3]);
	if (edgeTable[CubeIndex] & 8)
		VertexList[3] =
		VertexInterp(Grid.p[3],Grid.p[0],Grid.val[3],Grid.val[0]);
	if (edgeTable[CubeIndex] & 16)
		VertexList[4] =
		VertexInterp(Grid.p[4],Grid.p[5],Grid.val[4],Grid.val[5]);
	if (edgeTable[CubeIndex] & 32)
		VertexList[5] =
		VertexInterp(Grid.p[5],Grid.p[6],Grid.val[5],Grid.val[6]);
	if (edgeTable[CubeIndex] & 64)
		VertexList[6] =
		VertexInterp(Grid.p[6],Grid.p[7],Grid.val[6],Grid.val[7]);
	if (edgeTable[CubeIndex] & 128)
		VertexList[7] =
		VertexInterp(Grid.p[7],Grid.p[4],Grid.val[7],Grid.val[4]);
	if (edgeTable[CubeIndex] & 256)
		VertexList[8] =
		VertexInterp(Grid.p[0],Grid.p[4],Grid.val[0],Grid.val[4]);
	if (edgeTable[CubeIndex] & 512)
		VertexList[9] =
		VertexInterp(Grid.p[1],Grid.p[5],Grid.val[1],Grid.val[5]);
	if (edgeTable[CubeIndex] & 1024)
		VertexList[10] =
		VertexInterp(Grid.p[2],Grid.p[6],Grid.val[2],Grid.val[6]);
	if (edgeTable[CubeIndex] & 2048)
		VertexList[11] =
		VertexInterp(Grid.p[3],Grid.p[7],Grid.val[3],Grid.val[7]);

	//int NewVertexCount=0;
	for (int i=0;i<12;i++)
		LocalRemap[i] = -1;

	for (int i=0;triTable[CubeIndex][i]!=-1;i++)
	{
		if(LocalRemap[triTable[CubeIndex][i]] == -1)
		{
			NewVertexList[NewVertexCount] = VertexList[triTable[CubeIndex][i]];
			LocalRemap[triTable[CubeIndex][i]] = NewVertexCount;
			NewVertexCount++;
		}
	}

	vertices.resize(NewVertexCount);
	for (int i=0;i<NewVertexCount;i++) {
		vertices[i] = NewVertexList[i];
	}

	//cout<<LocalRemap[triTable[CubeIndex][0]]<<endl;
	TriangleCount = 0;
	for (int i=0;triTable[CubeIndex][i]!=-1;i+=3) {
		Vec3 a,b,c,a_norm,b_norm,c_norm;
		a = NewVertexList[LocalRemap[triTable[CubeIndex][i+0]]];
		b = NewVertexList[LocalRemap[triTable[CubeIndex][i+1]]];
		c = NewVertexList[LocalRemap[triTable[CubeIndex][i+2]]];
		//a_norm = normal_at_point(a);
		//b_norm = normal_at_point(b);
		//c_norm = normal_at_point(c);
		//        Vec3 midpoint = ((b-a) + (c-a))*.5f;
		//        Vec3 mid_normal = normal_at_point(midpoint);

		Triangle* new_tri = new Triangle(a,b,c);//a_norm,b_norm,c_norm);
		new_triangles.push_back(new_tri);
		TriangleCount++;
	}
	return new_triangles;
}

/******************
* Kernel Function *
******************/
//default is used for all terms except pressure and viscosity
float default_kernel(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	return (315/(64*PI*pow(H,9.0f)))*pow((H*H - mag),3.0f);
}

Vec3 default_gradient(Vec3 r_i,Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float coeff = -945.0f/(32.0f*PI*pow(H,9.0f));

	return diff_vec * coeff * pow((H*H - mag),2.0f);
}

float default_laplacian(Vec3 r_i,Vec3 r_j){

	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float coeff = -945.0f/(32.0f*PI*pow(H,9.0f));

	return coeff * (H*H - mag) * (3.0f*H*H - 7.0f*mag);
}

/*
Kernel used specifically for pressure gradient calculation.
*/
Vec3 pressure_kernel_gradient(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	float coeff = (45/(PI*pow(H,6.0f)))*pow((H-sqrt(mag)),2.0f);

	Vec3 v(0,0,0);
	//might want to return random vector if sqrt(mag)==0
	if(mag==0){
		cout<<'h'<<endl;
	}
	return diff_vec*(-coeff)/(sqrt(mag));
}

/*
Kernel used specifically for viscosity calculation.
*/
float viscosity_kernel_laplacian(Vec3 r_i, Vec3 r_j){
	Vec3 diff_vec = r_i-r_j;
	float mag = dot(diff_vec,diff_vec);

	return (45.0f/(PI*pow(H,6)))*(H-sqrt(mag));
}

/*******************
* Particle Methods *
*******************/
/*
Returns the density, using kernels, at the point param.
*/
float density_at_point(Vec3 point){
	//first should generate a list of the particles we need to check, using the box for this point.
	int box_number = NEIGHBOR.compute_box_num(point,SUPPORT_RADIUS,CONTAINER.max.x,CONTAINER.min.x);

	Particle* temp_particle;
	vector<int> neighbor_vec = NEIGHBOR.box_neighbors[box_number];
	float density = 0;//default_kernel(point,point);
	for (int i = 0; i<neighbor_vec.size(); i++){
		//update density. There will be a problem if the point is exactly equal to sum particle (divide by zero error).
		temp_particle = PARTICLES[neighbor_vec[i]];
		Vec3 diff_vec = point-temp_particle->position;
		float mag = dot(diff_vec,diff_vec);
		if(mag<H*H){
			density += temp_particle->mass*default_kernel(point,temp_particle->position);
		}
	}

	return density;
}

/*
Similarly returns the normal at the point param.
*/
Vec3 normal_at_point(Vec3 point){
	//set the normal at each point
	int box_number = NEIGHBOR.compute_box_num(point,SUPPORT_RADIUS,CONTAINER.max.x,CONTAINER.min.x);
	Particle* temp_particle;
	vector<int> neighbor_vec = NEIGHBOR.box_neighbors[box_number];
	Vec3 normal(0,0,0);
	for (int j = 0; j<neighbor_vec.size(); j++){
		temp_particle = PARTICLES[neighbor_vec[j]];
		Vec3 r = point-temp_particle->position;
		float mag = dot(r,r);
		if(mag<H*H){
			normal += default_gradient(point,temp_particle->position)*(temp_particle->mass)/(temp_particle->density);
		}

	}

	float length = sqrt(dot(normal,normal));
	if(length>0){
		normal = normal/length;
	}
	return normal*(-1.0f);//normal for triangles
}

/*
Create a new list of particles from the old list. Then, throw the old list.
To do this, calculate all quanities in Navier-Stokes, then use timestep to
update particle location from old location and velocity.
*/
void run_time_step(){
	vector<Particle*> new_particles; new_particles.resize(NUM_PARTICLES);
	vector<float> pressure_list; pressure_list.resize(NUM_PARTICLES);
	vector<Vec3> pressure_grad_list; pressure_grad_list.resize(NUM_PARTICLES);
	vector<Vec3> viscosity_list; viscosity_list.resize(NUM_PARTICLES);
	vector<float> color_list; color_list.resize(NUM_PARTICLES);
	vector<Vec3> tension_list; tension_list.resize(NUM_PARTICLES);

	Particle *base_particle, *temp_particle, *new_particle;
	float density = 0;

	//Sets density at each point
	for (int i = 0; i<NUM_PARTICLES; i++){
		density = 0;
		base_particle = PARTICLES[i];

		vector<int> neighbor_vec = base_particle->neighbors;
		int n = neighbor_vec.size();
		for (int j = 0; j<n; j++){ // changed to neighbors
			temp_particle = PARTICLES[neighbor_vec[j]];
			Vec3 r = base_particle->position-temp_particle->position;
			float mag = dot(r,r);
			if(mag<H*H){
				density += temp_particle->mass*default_kernel(base_particle->position,temp_particle->position);
			}
		}
		base_particle->density = density;
		pressure_list[i] = (STIFFNESS*(density-IDEAL_DENSITY));
	}

	//Sets pressure,viscosity laplacian, normal, and tension.
	for (int i = 0; i<NUM_PARTICLES; i++){
		base_particle = PARTICLES[i];

		Vec3 pressure_gradient(0,0,0);
		Vec3 viscosity_laplacian(0,0,0);
		float color = 0.0f;
		Vec3 normal(0,0,0);

		vector<int> neighbor_vec = base_particle->neighbors;
		for (int j = 0; j<neighbor_vec.size(); j++){ // changed to neighbors
			temp_particle = PARTICLES[neighbor_vec[j]];
			Vec3 r = base_particle->position-temp_particle->position;
			float mag = dot(r,r);
			if(mag<H*H){
				float weight = viscosity_kernel_laplacian(base_particle->position,temp_particle->position);
				viscosity_laplacian += ((temp_particle->velocity - base_particle->velocity)/temp_particle->density)*weight * temp_particle->mass;

				color += (temp_particle->mass / temp_particle->density) * default_laplacian(base_particle->position,temp_particle->position);

				normal += default_gradient(base_particle->position,temp_particle->position)*(temp_particle->mass / temp_particle->density);
				if(i!=neighbor_vec[j]){
					Vec3 weight_vec = pressure_kernel_gradient(base_particle->position,temp_particle->position);
					pressure_gradient += weight_vec * temp_particle->mass * ((pressure_list[i]+pressure_list[neighbor_vec[j]])/(2.0f*temp_particle->density)); 
				}
			}
		}
		pressure_grad_list[i] = (pressure_gradient*(-1.0f));
		viscosity_list[i] = (viscosity_laplacian);
		color_list[i] = (color);

		float length = sqrt(dot(normal,normal));
		if(length>TENSION_THRESHOLD){
			normal = normal/length;
			tension_list[i] = (normal*(-1.0f*color_list[i]));
		}else{
			Vec3 v(0,0,0);
			tension_list[i] = (v);
		}
	}

	//Create new particles from old particles and from pressure gradient.
	for (int i = 0; i<NUM_PARTICLES; i++){
		temp_particle = PARTICLES[i];

		//using Navier Stokes, calculate the change in velocity.

		Vec3 acceleration = GRAVITY + (tension_list[i]*SURFACE_TENSION
			+ (viscosity_list[i]*VISCOSITY)
			+ pressure_grad_list[i])/temp_particle->density;

		Vec3 position = temp_particle->position;
		Vec3 velocity = temp_particle->velocity;
		Vec3 new_velocity,new_position;

		new_velocity = velocity+acceleration*TIMESTEP;
		new_position = position + (velocity + acceleration*TIMESTEP*.5f)*TIMESTEP;

		float mass = temp_particle->mass;

		temp_particle = new Particle(new_position,new_velocity,mass,temp_particle->density);
		temp_particle->prev_position = position;

		CONTAINER.in_container(temp_particle,TIMESTEP); //applies reflections if outside of boundary.

		new_particles[i] = temp_particle;
	}

	PARTICLES = new_particles;

	//reset neighbor structure
    if (RENDERING_TRIANGLES) {
        // fill vector of vector of ints to access all particles that neighbor
        // a specific box number
        NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER,NUM_PARTICLES,true);
    } else {
        NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER,NUM_PARTICLES,false);
    }
}

/*
Code modified from Matthew Fisher. It is located at http://graphics.stanford.edu/~mdfisher/MarchingCubes.html
We have cited it in our final project write up
*/
void generate_triangles(){
	TRIANGLES.clear();
	GRIDCELL grid;

	//call polygonize, add triangles to TRIANGLES list.
	vector<Triangle*> new_triangles;
	for (float i = CONTAINER.min.x; i<CONTAINER.max.x; i = i+CUBE_TOL){
		for (float j = CONTAINER.min.y; j<CONTAINER.max.y; j = j+CUBE_TOL){
			for (float k = CONTAINER.min.z; k<CONTAINER.max.z; k = k+CUBE_TOL){
				new_triangles.clear();
				GRIDCELL grid;
				vector<Vec3> vertices;
				int vertex_count = 0;

				Vec3* vec_0 = new Vec3(i,j,k);
				grid.p.push_back(vec_0);
				grid.val.push_back(density_at_point(*vec_0)-DENSITY_TOL);

				Vec3* vec_1 = new Vec3(i+CUBE_TOL,j,k);
				grid.p.push_back(vec_1);
				grid.val.push_back(density_at_point(*vec_1)-DENSITY_TOL);

				Vec3* vec_2 = new Vec3(i+CUBE_TOL,j,k+CUBE_TOL);
				grid.p.push_back(vec_2);
				grid.val.push_back(density_at_point(*vec_2)-DENSITY_TOL);

				Vec3* vec_3 = new Vec3(i,j,k+CUBE_TOL);
				grid.p.push_back(vec_3);
				grid.val.push_back(density_at_point(*vec_3)-DENSITY_TOL);

				Vec3* vec_4 = new Vec3(i,j+CUBE_TOL,k);
				grid.p.push_back(vec_4);
				grid.val.push_back(density_at_point(*vec_4)-DENSITY_TOL);

				Vec3* vec_5 = new Vec3(i+CUBE_TOL,j+CUBE_TOL,k);
				grid.p.push_back(vec_5);
				grid.val.push_back(density_at_point(*vec_5)-DENSITY_TOL);

				Vec3* vec_6 = new Vec3(i+CUBE_TOL,j+CUBE_TOL,k+CUBE_TOL);
				grid.p.push_back(vec_6);
				grid.val.push_back(density_at_point(*vec_6)-DENSITY_TOL);

				Vec3* vec_7 = new Vec3(i,j+CUBE_TOL,k+CUBE_TOL);
				grid.p.push_back(vec_7);
				grid.val.push_back(density_at_point(*vec_7)-DENSITY_TOL);

				//make vertices and GRID object for marching cubes

				//call polygonize to make triangles
				new_triangles = polygonise(grid,vertex_count,vertices);

				//add triangles to full list
				TRIANGLES.insert(TRIANGLES.end(),new_triangles.begin(),new_triangles.end());
			}
		}
	}
}

/*****************
* OpenGL Methods *
*****************/
void initScene(){
	//create a list of random particles
	Particle* new_part;
	float noise = float(rand())/(float(RAND_MAX))*.1;
	float x,y,z,v_x,v_y,v_z;

    float step;
    // Sets up particles to preprogrammed scenes. Scene number taken from
    // command line options, or defaults to a scene.
    switch (SETUP_SCENE) {
        case 1:
            cout<<"Dam Break Scene"<<endl;
            step = .025;
            for(float i = 0; i<3.0f*(CONTAINER.max.x)/5.0f; i=i+step){
                for(float j = 0; j<2.0f*(CONTAINER.max.y)/5.0f; j=j+step){
                    for(float k = 0; k<3.0f*(CONTAINER.max.z)/5.0f; k=k+step){
                        Vec3 pos(i,j,k);
                        Vec3 vel(0,0,0);
                        new_part = new Particle(pos,vel,MASS,1000.0f);
                        PARTICLES.push_back(new_part);
                    }
                }
            }
            break;
            
        case 2:
            cout<<"Drop Scene"<<endl;
            step = .015;
            // 1.75, 3.25
            // 3, 4.5
            // 2, 3.5
            for(float i = 1.75*CONTAINER.max.x/5.0; i<3.25*(CONTAINER.max.x)/5.0f; i=i+step){
                for(float j = 3.0*CONTAINER.max.y/5.0f; j<4.5f*(CONTAINER.max.y)/5.0f; j=j+step){
                    for(float k = 2.0*CONTAINER.max.y/5.0f; k<4.0f*(CONTAINER.max.z)/5.0f; k=k+step){
                        Vec3 pos(i,j,k);
                        Vec3 vel(0,0,0);
                        new_part = new Particle(pos,vel,MASS,1000.0f);
                        PARTICLES.push_back(new_part);
                    }
                }
            }
            // some particles on floor
            for(float i = CONTAINER.min.x; i<(CONTAINER.max.x); i=i+step){
                for(float j = CONTAINER.min.y; j<.75f*(CONTAINER.max.y)/5.0f; j=j+step){
                    for(float k = CONTAINER.min.z; k<CONTAINER.max.z; k=k+step) {
                        noise = float(rand())/(float(RAND_MAX))*.05f;
                        Vec3 pos(i,j,k);
                        Vec3 vel(0,0,0);
                        new_part = new Particle(pos,vel,MASS,1000.f);
                        PARTICLES.push_back(new_part);
                    }
                }
            }
            break;
            
        case 3:
            cout<<"Velocity Towards Wall"<<endl;
            step = .025;
            for(float i = 0.1f*CONTAINER.max.x/5.0; i<5*(CONTAINER.max.x)/5.0f; i=i+step){
                for(float j = 0.1f*CONTAINER.max.y/5.0f; j<.5f*(CONTAINER.max.y)/5.0f; j=j+step){
                    for(float k = 0.1f*CONTAINER.max.y/5.0f; k<5.0f*(CONTAINER.max.z)/5.0f; k=k+step){
                        Vec3 pos(i,j*noise,k);
                        Vec3 vel(1.5f,0,0);
                        new_part = new Particle(pos,vel,MASS,1000.0f);
                        PARTICLES.push_back(new_part);
                    }
                }
            }
            break;
            
        case 4:
            cout<<"No Gravity Scene"<<endl;
            GRAVITY = Vec3(0,0,0);
            step = .02;
            for(float i = 2.0*CONTAINER.max.x/5.0; i<3.0f*(CONTAINER.max.x)/5.0f; i=i+step){
                for(float j = 2.5*CONTAINER.max.y/5.0f; j<3.5f*(CONTAINER.max.y)/5.0f; j=j+step){
                    for(float k = 1.5*CONTAINER.max.y/5.0f; k<4.0f*(CONTAINER.max.z)/5.0f; k=k+step){
                        //                        noise = float(rand())/(float(RAND_MAX))*.05f;
                        Vec3 pos(i,j,k);
                        Vec3 vel(0,0,0);
                        new_part = new Particle(pos,vel,MASS,1000.0f);
                        PARTICLES.push_back(new_part);
                    }
                }
            }
            break;
            
        case 5:
            cout<<"2 Glob Collision Scene"<<endl;
            step = .025;
            for(float i = 3.0*CONTAINER.max.x/5.0f; i<(CONTAINER.max.x); i=i+step){
                for(float j = 2.5*CONTAINER.max.y/4.0f; j<(CONTAINER.max.y); j=j+step){
                    for(float k = 2.0*CONTAINER.max.z/4.0f; k<(3.0f*CONTAINER.max.z/4.0f); k=k+step) {
                        noise = float(rand())/(float(RAND_MAX))*.05f;
                        Vec3 pos(i,j,k);
                        Vec3 vel(-1,-8,0);
                        new_part = new Particle(pos,vel,MASS,1000.0f);
                        PARTICLES.push_back(new_part);
                    }
                }
            }
            for(float i = CONTAINER.min.x; i<1.5f*(CONTAINER.max.x)/5.0f; i=i+step){
                for(float j = 2.5*CONTAINER.max.y/4.0f; j<(CONTAINER.max.y); j=j+step){
                    for(float k = 1.5*CONTAINER.max.z/4.0f; k<(3.0f*CONTAINER.max.z/4.0f); k=k+step) {
                        noise = float(rand())/(float(RAND_MAX))*.05f;
                        Vec3 pos(i,j,k);
                        Vec3 vel(5,-5,0);
                        new_part = new Particle(pos,vel,MASS,1000.f);
                        PARTICLES.push_back(new_part);
                    }
                }
            }
            // some particles on floor
            for(float i = CONTAINER.min.x; i<(CONTAINER.max.x); i=i+step){
                for(float j = CONTAINER.min.y; j<.5f*(CONTAINER.max.y)/5.0f; j=j+step){
                    for(float k = CONTAINER.min.z; k<CONTAINER.max.z; k=k+step) {
                        noise = float(rand())/(float(RAND_MAX))*.05f;
                        Vec3 pos(i,j,k);
                        Vec3 vel(0,0,0);
                        new_part = new Particle(pos,vel,MASS,1000.f);
                        PARTICLES.push_back(new_part);
                    }
                }
            }
            break;
            
        case 6:
            cout<<"Test Scene"<<endl;
            if (NUM_PARTICLES == 0) {
                // single particle test case
                Vec3 vel(0,0,0);
                Vec3 pos(.25,.25,.25);
                new_part = new Particle(pos,vel,MASS,1000);
                PARTICLES.push_back(new_part);
            }
            break;
            
        default:
            ////3D Uniform Scene
            cout<<"Default Horizontal Velocity Scene"<<endl;
            step = .025;
            for(float i = CONTAINER.min.x; i<(CONTAINER.max.x); i=i+step){
                for(float j = CONTAINER.min.y; j<1.0f*(CONTAINER.max.y)/5.0f; j=j+step){
                    for(float k = 1.0*CONTAINER.min.z; k<(CONTAINER.max.z); k=k+step){
//                        noise = float(rand())/(float(RAND_MAX))*.05f;
                        Vec3 pos(i,j,k);
                        Vec3 vel(-5,-3,0);
                        new_part = new Particle(pos,vel,MASS,1000.0f);
                        PARTICLES.push_back(new_part);
                    }
                }
            }
            break;
    }

	NUM_PARTICLES = PARTICLES.size();
	cout<<NUM_PARTICLES<<" Particles"<<endl;

    if (RENDERING_TRIANGLES) {
        NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER,NUM_PARTICLES,true);
    } else {
        NEIGHBOR.place_particles(PARTICLES,SUPPORT_RADIUS,CONTAINER,NUM_PARTICLES,false);
    }

	// create some lights
	GLfloat light_position[] = {1,1,1,0};
	GLfloat mat_specular[] = {0,0,0,1.0};
	GLfloat mat_diffuse[] = {0.0,0.0,1.0,1.0};
	GLfloat mat_ambient[] = {.1,.1,.1,1};
	GLfloat mat_shininess[] = {20.0};

    // set material properties
	glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
	glMaterialfv(GL_FRONT,GL_SHININESS,mat_shininess);
	glMaterialfv(GL_FRONT,GL_DIFFUSE,mat_diffuse);
	glMaterialfv(GL_FRONT,GL_AMBIENT,mat_ambient);
	glLightfv(GL_LIGHT0,GL_POSITION,light_position);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glClearColor(1,1,1,1);
}

void myDisplay(){
    if (HIDE) { glutHideWindow(); }
    
	float time_start = glutGet(GLUT_ELAPSED_TIME);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0,0,PIC_WIDTH,PIC_HEIGHT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(50,1.0f,.0001,1000);
	//glOrtho(CONTAINER.min.x,CONTAINER.max.x,CONTAINER.min.y,CONTAINER.max.y,CONTAINER.min.z,CONTAINER.max.z);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    // view for 0 to .5 cube
	//gluLookAt(.25f,.4f,1.2f,.25f,.18f,0.0f,0,1,0);
    
    // view for 0 to 2 cube
    gluLookAt(1.f,1.8f,4.5f, 1.f,.4f,0, 0,1,0);

    if (!PAUSE) {
        run_time_step();
        CURRENT_TIME += TIMESTEP;
    }
    
	if(CURRENT_TIME>LIFETIME){
		exit(0);
	}

	glEnable(GL_LIGHTING);
	if(RENDERING_TRIANGLES){
		glEnable(GL_LIGHTING);
		////draw triangles
		generate_triangles();
		Triangle *temp_triangle;
		for (int i = 0; i<TRIANGLES.size(); i++){
			temp_triangle = TRIANGLES[i];

			glClearColor(1,1,1,1);
			//glColor3f(1.0f,1.0f,1.0f);

			// wireframe
            if (RENDERING_WIREFRAME) {
                
                /*renders object twice as white and black with an offset. The offset removes some
                 artifacts from lines behind object.
                 */
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                
                glDisable(GL_LIGHTING);
                glClearColor(1,1,1,1);
                glColor3f(.4f,.4f,.4f);//color as white for wireframe
                
                glBegin(GL_TRIANGLES);
                glVertex3f(temp_triangle->a.x,temp_triangle->a.y,temp_triangle->a.z);
                glNormal3f(temp_triangle->a_normal.x,temp_triangle->a_normal.y,temp_triangle->a_normal.z);
                glVertex3f(temp_triangle->b.x,temp_triangle->b.y,temp_triangle->b.z);
                glNormal3f(temp_triangle->b_normal.x,temp_triangle->b_normal.y,temp_triangle->b_normal.z);
                glVertex3f(temp_triangle->c.x,temp_triangle->c.y,temp_triangle->c.z);
                glNormal3f(temp_triangle->c_normal.x,temp_triangle->c_normal.y,temp_triangle->c_normal.z);
                glEnd();
                
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//color as background color.
                glEnable(GL_POLYGON_OFFSET_FILL);
                glPolygonOffset(1.0,1.0);
                glClearColor(1,1,1,1);
                glColor3f(1,1,1);
                
                glBegin(GL_TRIANGLES);
                glVertex3f(temp_triangle->a.x,temp_triangle->a.y,temp_triangle->a.z);
                glNormal3f(temp_triangle->a_normal.x,temp_triangle->a_normal.y,temp_triangle->a_normal.z);
                glVertex3f(temp_triangle->b.x,temp_triangle->b.y,temp_triangle->b.z);
                glNormal3f(temp_triangle->b_normal.x,temp_triangle->b_normal.y,temp_triangle->b_normal.z);
                glVertex3f(temp_triangle->c.x,temp_triangle->c.y,temp_triangle->c.z);
                glNormal3f(temp_triangle->c_normal.x,temp_triangle->c_normal.y,temp_triangle->c_normal.z);
                glEnd();
                
                glDisable(GL_POLYGON_OFFSET_FILL);
                
            } else {
                glBegin(GL_TRIANGLES);
                glVertex3f(temp_triangle->a.x,temp_triangle->a.y,temp_triangle->a.z);
                glNormal3f(temp_triangle->a_normal.x,temp_triangle->a_normal.y,temp_triangle->a_normal.z);
                glVertex3f(temp_triangle->b.x,temp_triangle->b.y,temp_triangle->b.z);
                glNormal3f(temp_triangle->b_normal.x,temp_triangle->b_normal.y,temp_triangle->b_normal.z);
                glVertex3f(temp_triangle->c.x,temp_triangle->c.y,temp_triangle->c.z);
                glNormal3f(temp_triangle->c_normal.x,temp_triangle->c_normal.y,temp_triangle->c_normal.z);
                glEnd();
            }
		}
	}else{
		//draw particles
		glPointSize(4.0f);
		//glBegin(GL_POINTS);
		glEnable(GL_LIGHTING);
		Particle* temp_part;
		glClearColor(1,1,1,1);
		for (int i = 0; i<PARTICLES.size(); i++){
			temp_part = PARTICLES[i];

			//Draw sphere of radius H around particles
			//glColor3f(0,0,1.0);

			//glVertex3f(temp_part->position.x,temp_part->position.y,temp_part->position.z);
			glPushMatrix();
			glTranslated(temp_part->position.x,temp_part->position.y,temp_part->position.z);
			glutSolidSphere(DRAW_RADIUS,16,16);
			glPopMatrix();

		}
		//glEnd();
	}

	//Draw floor
	glDisable(GL_LIGHTING);
	glClearColor(1,1,1,1);
	glColor3f(.7f,.7f,.7f);
	glBegin(GL_POLYGON);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glNormal3f(0,1,0);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glNormal3f(0,1,0);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glNormal3f(0,1,0);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glNormal3f(0,1,0);
	glEnd();

	//Draw wireframe container
	glPolygonMode(GL_FRONT, GL_LINE);
	glPolygonMode(GL_BACK, GL_LINE);

	glDisable(GL_LIGHTING);
	glClearColor (1, 1, 1, 1);
	glColor3f(.8f,.8f,.8f);

	glBegin(GL_POLYGON);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.max.y+DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.max.y+DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.max.y+DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.max.y+DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.max.y+DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glVertex3f(CONTAINER.min.x-DRAW_RADIUS,CONTAINER.max.y+DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.max.y+DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.max.y+DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.min.z-DRAW_RADIUS);
	glVertex3f(CONTAINER.max.x+DRAW_RADIUS,CONTAINER.min.y-DRAW_RADIUS,CONTAINER.max.z+DRAW_RADIUS);
	glEnd();

	glPolygonMode(GL_FRONT, GL_FILL); // fill mode
	glPolygonMode(GL_BACK, GL_FILL);

	glPopMatrix();

	if ((OUTPUT_IMAGE && IMAGE_DELAY == 0 )|| OUTPUT_SINGLE_IMAGE) {
		// Output image to file
		FreeImage_Initialise();

		// Make the BYTE array, factor of 3 because it's RBG.
		BYTE* pixels = new BYTE[ 3 * PIC_HEIGHT * PIC_WIDTH];

		glReadPixels(0, 0, PIC_WIDTH, PIC_HEIGHT, GL_BGR, GL_UNSIGNED_BYTE, pixels);

		// Convert to FreeImage format & save to file
		FIBITMAP* image = FreeImage_ConvertFromRawBits(pixels, PIC_WIDTH, PIC_HEIGHT, 3 * PIC_WIDTH, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);

		std::stringstream ss;
		ss << IMAGE_COUNTER;
		std::string s(ss.str());
		string name = std::string("images/")+s+".png";
		FreeImage_Save(FIF_PNG, image, name.c_str(), 0);

		// FreeImage resources
		FreeImage_Unload(image);
		delete [] pixels;

        // Increment image counter, terminate if 350 frames reached
		IMAGE_COUNTER++;
		if (IMAGE_COUNTER == 350) {
			exit(0);
		}
		if (OUTPUT_SINGLE_IMAGE) {
			OUTPUT_SINGLE_IMAGE = false;
			cout<<"Saved single image."<<endl;
		}
	}

    // print out the frames per second to screen
	if (SHOW_FPS) {
		float time_end = glutGet(GLUT_ELAPSED_TIME);
		cout<<1000/(time_end - time_start)<<" FPS"<<endl;
		SHOW_FPS = false;
	}

    // Raytrace Movie Call
    if (RAYTRACE_MOVIE && IMAGE_DELAY == 0) {
        output_obj();
    }
    if (IMAGE_DELAY != 0) {
        IMAGE_DELAY--;
        IMAGE_COUNTER++;
    }
    
	glFlush();
	glutSwapBuffers();
}

int main(int argc, char* argv[]){

	// Command Line Arg Parsing
	for (int i = 1; i < argc; ) {
		if (strcmp(argv[i],"-m") == 0) {
			cout<<"Generating images for animation."<<endl;
			OUTPUT_IMAGE = true;
			i += 1;
			continue;
		}
		if (strcmp(argv[i],"-s") == 0) {
			SETUP_SCENE = atoi(argv[i+1]);
			i += 2;
			continue;
		}
		if (strcmp(argv[i],"-help") == 0) {
			cout<<"Help Info for UC Berkeley CS184 Spring 2013 Final Project"<<endl;
			cout<<"by Tyler Brabham and Zack Mayeda"<<endl;
			cout<<endl;
			cout<<"Command Line Arguments:"<<endl;
            cout<<"-help : List all options"<<endl;
            cout<<"-s # : Load scene number #"<<endl;
            cout<<"       1 - Dam Break"<<endl;
            cout<<"       2 - Drop"<<endl;
            cout<<"       3 - Randomized"<<endl;
            cout<<"       4 - No Gravity"<<endl;
            cout<<"       5 - Collision"<<endl;
            cout<<"       Default - Horizontal Velocity"<<endl;
            cout<<"-m : Export Still Frames for movie into Images/"<<endl;
            cout<<"-rm : Export Raytraced Still Frames for movie into Multi_Trace/output_pics/"<<endl;
            cout<<"-i : Render Isosurface"<<endl;
            cout<<"-p : Render Particles"<<endl;
            cout<<"-v : Increase viscosity and surface tension"<<endl;
            cout<<"-h : Hide Display Window"<<endl;
            cout<<"-delay # : delay movie output by # frames"<<endl;
            cout<<"======================="<<endl;
			cout<<"Live Commands:"<<endl;
            cout<<"'i' : toggle isosurface/particles"<<endl;
            cout<<"'r' : output single raytraced image and exit"<<endl;
            cout<<"'f' : output current FPS to command line"<<endl;
            cout<<"'p' : save current image, not raytrace, and continue"<<endl;
            cout<<"spacebar : pause"<<endl;
            cout<<"'q' : quit"<<endl;
			i += 1;
			exit(0);
		}
        if (strcmp(argv[i],"-rm") == 0) {
            cout<<"Outputting raytrace OBJ files for movie."<<endl;
            RAYTRACE_MOVIE = true;
            i += 1;
            continue;
        }
        if (strcmp(argv[i],"-i") == 0) {
            RENDERING_TRIANGLES = true;
            i += 1;
            continue;
        }
        if (strcmp(argv[i],"-p") == 0) {
            RENDERING_TRIANGLES = false;
            i += 1;
            continue;
        }
        if (strcmp(argv[i],"-delay") == 0) {
            IMAGE_DELAY = atoi(argv[i+1]);
            i += 2;
            continue;
        }
        if (strcmp(argv[i],"-w") == 0) {
            RENDERING_WIREFRAME = true;
            CUBE_TOL = .025f;
            i += 1;
            continue;
        }
        if (strcmp(argv[i],"-v") == 0) {
            VISCOSITY = 30.f;
            SURFACE_TENSION = 10.0f;
            i += 1;
            continue;
        }
        if (strcmp(argv[i],"-h") == 0) {
            HIDE = true;
            i += 1;
            continue;
        }
	}

	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

	//The size and position of the window
	glutInitWindowSize(PIC_WIDTH,PIC_HEIGHT);
	glutInitWindowPosition(0,0);
	glutCreateWindow("Tyler and Zack Final Project");
    //glutHideWindow();
    
	initScene();

	glutDisplayFunc(myDisplay);
	glutIdleFunc(myDisplay);
	glutKeyboardFunc(keyPressed);

	glutMainLoop();

	return 0;
}
