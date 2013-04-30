#include "../Raytracer/Raytracer.h"

#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]) {
    char* input_filename;
    char* output_filename;
    int height, width;
    
    if (argc > 5) {
        input_filename = argv[2];
        output_filename = argv[3];
        width = atoi(argv[4]);
        height = atoi(argv[5]);
    } else {
        cout<<"Error, not called with correct arguments."<<endl;
        cout<<"./start_trace input_file output_file width height"<<endl;
        return 1;
    }
    
    Raytracer r;
    r.ray_trace_start(input_filename,output_filename,width,height);
    
    return 0;
}