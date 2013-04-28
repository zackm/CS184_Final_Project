#pragma once
#include "MarchingCubes.h"

#include "Particle.h"

#include "Triangle.h"

#include <vector>

#include <GL/glut.h>
#include <GL/glu.h>

#include <stdlib.h>
#include <iostream>

#include <cmath>

#include "Vec3.h"

#include "Container.h"

#include "Neighbor.h"

#include "Raytracer/Raytracer.h"

struct GRIDCELL{
   vector<Vec3*> p;		//position of each corner of the grid in world space
   vector<float> val;	//value of the function at this grid corner
};