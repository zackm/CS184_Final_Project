#include "../Neighbor.h"

class Raytracer {
public:
    int ray_trace_start(std::string,std::string,int,int);
	Neighbor neighborhood;
    
	Raytracer(Neighbor neighbor_arg){neighborhood = neighbor_arg;};
	Raytracer(){};
};
