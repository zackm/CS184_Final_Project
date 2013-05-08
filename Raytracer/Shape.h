#pragma once
#include "Ray.h"

#pragma once
#include "LocalGeo.h"

#pragma once
#include "BRDF.h"

#include "Transformation.h"

using namespace std;

class Shape {
public:
	BRDF brdf;
	float index_of_refraction;
	bool transparency;

	Shape(void){};
	~Shape(void){};
	virtual bool intersect(Ray&, float*, LocalGeo*,bool*) =0;
	virtual bool intersect(Ray&,bool*) =0;
	virtual BRDF get_brdf() =0;
};

