/*
CS184 Assignment 2 - Ray Tracing
Tyler Brabham cs184-ej
Zack Mayeda cs184-bg
*/

#include "glm/glm.hpp"
// OpenGL Math Library
// http://glm.g-truc.net/code.html
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

#include <list>

// look into pragma once
#ifndef __FILM_H__
#define __FILM_H__
#include "Film.h"
#endif

#ifndef __CAMERA_H__
#define __CAMERA_H__
#include "Camera.h"
#endif

#include "Scene.h"

#include "Sphere.h"

#include "DirectionalLight.h"

#include "PointLight.h"

#include <vector>

#include "RTriangle.h"

#include "Transformation.h"

#include "Raytracer.h"

using namespace std;

const float pi = 3.14159265359;

bool OBJ_ON = false; // flag for OBJ input file parsing

bool flag = true; // temp flag for creating one sphere

/*
Simply creates the 4 by 4 rotation matrix where [x,y,z]
is the axis of rotation and angle is the angle in degrees
to rotate.
*/
glm::mat4 create_rotate(float x, float y, float z, float angle) {
	glm::mat4 rotate_mat;
	angle = angle * pi / 180.0f;
	glm::mat4 id(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
	glm::mat4 cos_mat(x*x,x*y,x*z,0,x*y,y*y,y*z,0,x*z,y*z,z*z,0,0,0,0,0);
	glm::mat4 sin_mat(0,z,-y,0,-z,0,x,0,y,-x,0,0,0,0,0,0);
	rotate_mat = glm::cos(angle) * id + (1 - glm::cos(angle)) * cos_mat + glm::sin(angle) * sin_mat;
	rotate_mat[3][3] = 1;
	return rotate_mat;
}

/*
Helper method for obj file parsing. Finds the number of slashes in
a face reference. This is due to multiple standards for .obj file format.
*/
int slash_counter(string s) {
	int count = 0;
	int pos = 0;
	while (s.find("/", pos) != std::string::npos) {
		count++;
		int curr_pos = s.find("/",pos);
		pos = curr_pos + 1;
	}
	return count;
}

/*
Method to set the camera position using the max and min points in the object.
*/
void set_camera_and_perspective(Camera c, glm::vec3 max, glm::vec3 min) {
	glm::vec3 center,camera_pos,camera_up;

	float diameter = glm::max(glm::max(max.x-min.x,max.z-min.z),max.y-min.y);

	center.x = (max.x+min.x)/2.0f;
	center.y = (max.y+min.y)/2.0f;
	center.z = (max.z+min.z)/2.0f;
	camera_pos.x = camera_pos.y = 0;
	camera_pos.z = 0;

	camera_up.x = camera_up.z = 0;
	camera_up.y = 1.0f;

	float fov = 90;

	center.x = 0; center.y = 0, center.z = -3;
	// from, to .up
	c.set_args(camera_pos,center,camera_up,fov);

	cout<<"Max:"<<max.x<<", "<<max.y<<", "<<max.z<<endl;
	cout<<"Min:"<<min.x<<", "<<min.y<<", "<<min.z<<endl;
	cout<<"camera pos = "<<camera_pos.x<<", "<<camera_pos.y<<", "<<camera_pos.z<<endl;
	cout<<"center = "<<center.x<<", "<<center.y<<", "<<center.z<<endl;
	// exit(0);
}

int Raytracer::ray_trace_start() {
	int WIDTH = 1200;
	int HEIGHT = 1200;

	Scene s;
	Camera c;
	int maxdepth = 5;
	std::string output_name;
	vector<glm::vec3> vertices;
	vector<glm::vec3> vertexnorm_v;
	vector<glm::vec3> vertexnorm_n;
	vector<glm::mat4> mat_stack;
	glm::mat4 current_mat;

	glm::vec3 max(0,0,0); // need to set these to max float, min float probably
	glm::vec3 min(0,0,0);

	//default material properties
	glm::vec3 ka(.2f, .2f, .2f);
	glm::vec3 kd(0,0,0);
	glm::vec3 ks(0,0,0);
	glm::vec3 kr(0,0,0);
	glm::vec3 ke(0,0,0);
	float sp = 1;

	bool has_reflect_coeff = false;

	// Initializers for AS2 File Format
	// Push identity matrix onto stack
	glm::mat4 id_mat(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
	mat_stack.push_back(id_mat); // top of the stack is the end of the list

	// Initializers for OBJ Format
	vector<glm::vec3> vert_list;
	vector<glm::vec3> norm_list;
	bool OBJ_NORM = false;

	// fill first position for proper numbering, since obj parsing starts at  instead of 0.
	vert_list.push_back(glm::vec3(0,0,0));
	norm_list.push_back(glm::vec3(0,0,0));
	vector<RTriangle*> filler_tri;

	// Filename Business
//	if (argc < 2) {
//		cout << "No filname given. Terminating" << endl;
//		exit(1);
//	}

	//std::string filename = argv[1];
    std::string filename = "fluid.obj";
	cout << "Filename " << filename << " found." << endl;
	if (filename.find(".obj") != std::string::npos) {
		cout<<"OBJ Input File Detected."<<endl;
		OBJ_ON = true;
		// Hard code a light for OBJ
		float x = 1;
		float y = 1;
		float z = 6;
		float r = .9;
		float g = .9;
		float b = .9;
		Transformation directional_trans(mat_stack);
		DirectionalLight* dl = new DirectionalLight(glm::vec3(x,y,z),glm::vec3(r,g,b),directional_trans);
		s.add_light(dl);

	}

	// Arg Parser
	std::ifstream inpfile(filename.c_str());
	if(!inpfile.is_open()) {
		std::cout << "Unable to open file" << std::endl;
	} else {
		std::string line;
		//MatrixStack mst;
		while(inpfile.good()) {
			std::vector<std::string> splitline;
			std::string buf;
			std::getline(inpfile,line);
			std::stringstream ss(line);
			while (ss >> buf) {
				splitline.push_back(buf);
			}
			//Ignore blank lines
			if(splitline.size() == 0) {
				continue;
			}
			//Ignore comments
			if(splitline[0][0] == '#') {
				continue;
			}

			// Non-OBJ input files (AS2 format)
			if (!OBJ_ON) {
				//Valid commands:
				//size width height
				//  must be first command of file, controls image size
				if(!splitline[0].compare("size")) {
					WIDTH = atoi(splitline[1].c_str());
					HEIGHT = atoi(splitline[2].c_str());
				}

				//maxdepth depth
				//  max # of bounces for ray (default 5)
				else if(!splitline[0].compare("maxdepth")) {
					maxdepth = atoi(splitline[1].c_str());
				}

				//output filename
				//  output file to write image to 
				else if(!splitline[0].compare("output")) {
					output_name = splitline[1];
				}	

				//sphere x y z radius
				//  Deﬁnes a sphere with a given position and radius.
				else if(!splitline[0].compare("sphere")) {
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					float r = atof(splitline[4].c_str());
					// Create new sphere:
					//   Store 4 numbers
					//   Store current property values
					//   Store current top of matrix stack

					//make transformation matrix

					Transformation sphere_trans(mat_stack);
					Sphere* sph = new Sphere(glm::vec3(x,y,z),r,ka,kd,ks,kr,ke,sp,sphere_trans);
					s.add_shape(sph);
				}

				//vertex x y z
				//  Deﬁnes a vertex at the given location.
				//  The vertex is put into a pile, starting to be numbered at 0.
				else if(!splitline[0].compare("vertex")) {
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					// Create a new vertex with these 3 values, store in some array
					glm::vec3 vert(x,y,z);
					vertices.push_back(vert);
				}

				//vertexnormal x y z nx ny nz
				//  Similar to the above, but deﬁne a surface normal with each vertex.
				//  The vertex and vertexnormal set of vertices are completely independent
				//  (as are maxverts and maxvertnorms).
				else if(!splitline[0].compare("vertexnormal")) {
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					float nx = atof(splitline[4].c_str());
					float ny = atof(splitline[5].c_str());
					float nz = atof(splitline[6].c_str());
					// Create a new vertex+normal with these 6 values, store in some array
					glm::vec3 norm_v(x,y,z);
					glm::vec3 norm_n(nx,ny,nz);
					vertexnorm_v.push_back(norm_v);
					vertexnorm_n.push_back(norm_n);
				}

				//tri v1 v2 v3
				//  Create a triangle out of the vertices involved (which have previously been speciﬁed with
				//  the vertex command). The vertices are assumed to be speciﬁed in counter-clockwise order. Your code
				//  should internally compute a face normal for this triangle.
				else if(!splitline[0].compare("tri")) {
					int v1 = atoi(splitline[1].c_str());
					int v2 = atoi(splitline[2].c_str());
					int v3 = atoi(splitline[3].c_str());
					// Create new triangle:
					//   Store pointer to array of vertices
					//   Store 3 integers to index into array
					//   Store current property values
					//   Store current top of matrix stack

					Transformation tri_trans(mat_stack);
					glm::vec3 vert_1 = tri_trans.world_point(vertices[v1]);
					glm::vec3 vert_2 = tri_trans.world_point(vertices[v2]);
					glm::vec3 vert_3 = tri_trans.world_point(vertices[v3]);

					//RTriangle *t = new RTriangle(vert_1,vert_2,vert_3,ka,kd,ks,kr,ke,sp);
					//t->index_of_refraction = 1.33f;
					//t->transparency = true;//true;//not really working for triangles
					//s.add_shape(t);
				}

				//trinormal v1 v2 v3
				//  Same as above but for vertices speciﬁed with normals.
				//  In this case, each vertex has an associated normal, 
				//  and when doing shading, you should interpolate the normals 
				//  for intermediate points on the triangle.
				else if(!splitline[0].compare("trinormal")) {
					int v1 = atoi(splitline[1].c_str());
					int v2 = atoi(splitline[2].c_str());
					int v3 = atoi(splitline[3].c_str());
					// Create new triangle:
					//   Store pointer to array of vertices (Different array than above)
					//   Store 3 integers to index into array
					//   Store current property values
					//   Store current top of matrix stack

					Transformation tri_trans(mat_stack);
					glm::vec3 vert_1 = tri_trans.world_point(vertexnorm_v[v1]);
					glm::vec3 vert_2 = tri_trans.world_point(vertexnorm_v[v2]);
					glm::vec3 vert_3 = tri_trans.world_point(vertexnorm_v[v3]);

					glm::vec3 norm_1 = tri_trans.world_normal(vertexnorm_n[v1]);
					glm::vec3 norm_2 = tri_trans.world_normal(vertexnorm_n[v2]);
					glm::vec3 norm_3 = tri_trans.world_normal(vertexnorm_n[v3]);

					RTriangle *t = new RTriangle(vert_1,vert_2,vert_3,ka,kd,ks,kr,ke,sp,
						norm_1,norm_2,norm_3);
					s.add_shape(t);
				}

				//directional x y z r g b
				//  The direction to the light source, and the color, as in OpenGL.
				else if(!splitline[0].compare("directional")) {
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					float r = atof(splitline[4].c_str());
					float g = atof(splitline[5].c_str());
					float b = atof(splitline[6].c_str());

					Transformation directional_trans(mat_stack);
					DirectionalLight* dl = new DirectionalLight(glm::vec3(x,y,z),glm::vec3(r,g,b),directional_trans);
					s.add_light(dl);
				}

				//point x y z r g 
				//  The location of a point source and the color, as in OpenGL.
				else if(!splitline[0].compare("point")) {
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					float r = atof(splitline[4].c_str());
					float g = atof(splitline[5].c_str());
					float b = atof(splitline[6].c_str());

					Transformation point_trans(mat_stack);
					PointLight* pt = new PointLight(glm::vec3(x,y,z),glm::vec3(r,g,b),point_trans);
					s.add_light(pt);
				}

				// camera lookfromx lookfromy lookfromz lookatx lookaty lookatz upx upy upz fov
				else if(!splitline[0].compare("camera")){
					float from_x = atof(splitline[1].c_str());
					float from_y = atof(splitline[2].c_str());
					float from_z = atof(splitline[3].c_str());
					float to_x = atof(splitline[4].c_str());
					float to_y = atof(splitline[5].c_str());
					float to_z = atof(splitline[6].c_str());
					float up_x = atof(splitline[7].c_str());
					float up_y = atof(splitline[8].c_str());
					float up_z = atof(splitline[9].c_str());
					float fov = atof(splitline[10].c_str());

					Camera cam(glm::vec3(from_x,from_y,from_z),glm::vec3(to_x,to_y,to_z),glm::vec3(up_x,up_y,up_z),fov);
					c = cam;
				}

				//ambient r g b
				//  The global ambient color to be added for each object 
				//  (default is .2,.2,.2)
				else if(!splitline[0].compare("ambient")) {
					float r = atof(splitline[1].c_str());
					float g = atof(splitline[2].c_str());
					float b = atof(splitline[3].c_str());
					ka = glm::vec3(r,g,b);
				}

				//diﬀuse r g b
				//  speciﬁes the diﬀuse color of the surface.
				else if(!splitline[0].compare("diffuse")) {
					float r = atof(splitline[1].c_str());
					float g = atof(splitline[2].c_str());
					float b = atof(splitline[3].c_str());
					kd = glm::vec3(r,g,b);
				}

				//specular r g b 
				//  speciﬁes the specular color of the surface.
				else if(!splitline[0].compare("specular")) {
					float r = atof(splitline[1].c_str());
					float g = atof(splitline[2].c_str());
					float b = atof(splitline[3].c_str());
					ks = glm::vec3(r,g,b);

					if (!has_reflect_coeff){
						kr = ks;
					}
				}

				//reflection r g b 
				//  speciﬁes the reflective color of the surface.
				else if(!splitline[0].compare("reflect")) {
					float r = atof(splitline[1].c_str());
					float g = atof(splitline[2].c_str());
					float b = atof(splitline[3].c_str());
					kr = glm::vec3(r,g,b);
					has_reflect_coeff = true;
					// unsure of command line arg ??
				}			

				else if(!splitline[0].compare("emission")) {
					float r = atof(splitline[1].c_str());
					float g = atof(splitline[2].c_str());
					float b = atof(splitline[3].c_str());
					ke = glm::vec3(r,g,b);
				}

				//shininess s
				//  speciﬁes the shininess of the surface.
				else if(!splitline[0].compare("shininess")) {
					sp = atof(splitline[1].c_str());
				} 

				//translate x y z
				//  A translation 3-vector
				else if(!splitline[0].compare("translate")) {
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					// Update top of matrix stack
					// Matrix of form
					// 1 0 0 tx
					// 0 1 0 ty
					// 0 0 1 tz
					// 0 0 0 1

					glm::mat4 translate_mat(1,0,0,0,0,1,0,0,0,0,1,0,x,y,z,1);
					current_mat = mat_stack.back();

					mat_stack.pop_back();
					current_mat = current_mat * translate_mat;

					mat_stack.push_back(current_mat);
				}

				//rotate x y z angle
				//  Rotate by angle (in degrees) about the given axis as in OpenGL.
				else if(!splitline[0].compare("rotate")) {
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					float angle = atof(splitline[4].c_str());
					// Update top of matrix stack

					glm::mat4 rotate_mat = create_rotate(x,y,z,angle);
					current_mat = mat_stack.back();

					mat_stack.pop_back();
					current_mat = current_mat * rotate_mat;

					mat_stack.push_back(current_mat);
				}

				//scale x y z
				//  Scale by the corresponding amount in each axis (a non-uniform scaling).
				else if(!splitline[0].compare("scale")) {
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					// Update top of matrix stack
					// Matrix of form
					// sx 0  0  0
					// 0  sy 0  0
					// 0  0  sz 0
					// 0  0  0  1
					glm::mat4 scale_mat(x,0,0,0,0,y,0,0,0,0,z,0,0,0,0,1);
					current_mat = mat_stack.back();

					mat_stack.pop_back();
					current_mat = current_mat * scale_mat;

					mat_stack.push_back(current_mat);
				}

				//pushTransform
				//  Push the current modeling transform on the stack as in OpenGL. 
				//  You might want to do pushTransform immediately after setting 
				//   the camera to preserve the “identity” transformation.
				else if(!splitline[0].compare("pushTransform")) {
					//Transformation t_copy(current_trans.m);

					current_mat = glm::mat4(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
					mat_stack.push_back(current_mat);
				}

				//popTransform
				//  Pop the current transform from the stack as in OpenGL. 
				//  The sequence of popTransform and pushTransform can be used if 
				//  desired before every primitive to reset the transformation 
				//  (assuming the initial camera transformation is on the stack as 
				//  discussed above).
				else if(!splitline[0].compare("popTransform")) {
					mat_stack.pop_back();
				} else {
					std::cerr << "Unknown command: " << splitline[0] << std::endl;
				}
			} else {
				// Parse OBJ file
				if(!splitline[0].compare("s")){
					// ignore s commands
					continue;
				}
				else if(!splitline[0].compare("g")){
					// ignore g commands
					continue;
				}
				else if(!splitline[0].compare("vt")){
					// ignore vertex texture commands
					continue;
				}
				else if(!splitline[0].compare("usemtl")){
					// ignore material commands
					continue;
				}
				else if(!splitline[0].compare("v")){
					// add vertex to list
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					vert_list.push_back(glm::vec3(x,y,z));

					// compare current vertex with max and min
					if (x > max.x && y > max.y && z > max.z) {
						max.x = x;
						max.y = y;
						max.z = z;
					}
					if (x < min.x && y < min.y && z < min.z) {
						min.x = x;
						min.y = y;
						min.z = z;
					}
				}
				else if(!splitline[0].compare("vn")) {
					// add vertex normal to list
					OBJ_NORM = true;
					float x = atof(splitline[1].c_str());
					float y = atof(splitline[2].c_str());
					float z = atof(splitline[3].c_str());
					norm_list.push_back(glm::vec3(x,y,z));
				}
				//Multiple kinds of face formats. This handles a few of them.
				else if(!splitline[0].compare("f")) {
					int a_point,b_point,c_point,a_point_norm,b_point_norm,c_point_norm;
					glm::vec3 a,b,c;
					// add triangle to list, four different parsing options:
					// v for vertex, t for vextex texture (ignored), n for normal
					// face with vertecies v
					// face with vertex textures v/t
					// face with vertex norms v//n
					// face with text and norms v/t/n
					int count = slash_counter(splitline[1]);
					if (count == 1) {
						// case: f v/t
						OBJ_NORM = false;
						int pos = splitline[1].find("/");
						a_point = atoi(splitline[1].substr(0,pos).c_str());
						pos = splitline[2].find("/");
						b_point = atoi(splitline[2].substr(0,pos).c_str());
						pos = splitline[3].find("/");
						c_point = atoi(splitline[3].substr(0,pos).c_str());
					}
					else if (splitline[1].find("//") != std::string::npos) {
						// case: f v//n
						OBJ_NORM = true;
						int pos = splitline[1].find("/");
						a_point = atoi(splitline[1].substr(0,pos).c_str());
						a_point_norm = atoi(splitline[1].substr(pos+2).c_str());
						pos = splitline[2].find("/");
						b_point = atoi(splitline[2].substr(0,pos).c_str());
						b_point_norm = atoi(splitline[2].substr(pos+2).c_str());

						pos = splitline[3].find("/");
						c_point = atoi(splitline[3].substr(0,pos).c_str());
						c_point_norm = atoi(splitline[3].substr(pos+2).c_str());

					}
					else if (count == 2) {
						// case: f v/t/n
						OBJ_NORM = true;
						int pos = splitline[1].find("/");
						a_point = atoi(splitline[1].substr(0,pos).c_str());
						pos = splitline[1].find("/",pos+1);
						a_point_norm = atoi(splitline[1].substr(pos+1).c_str());

						pos = splitline[2].find("/");
						b_point = atoi(splitline[2].substr(0,pos).c_str());
						pos = splitline[2].find("/",pos+1);
						b_point_norm = atoi(splitline[2].substr(pos+1).c_str());

						pos = splitline[3].find("/");
						c_point = atoi(splitline[3].substr(0,pos).c_str());
						pos = splitline[3].find("/",pos+1);
						c_point_norm = atoi(splitline[3].substr(pos+1).c_str());

					} else {
						// case: f v
						OBJ_NORM = false;
						a_point = atoi(splitline[1].c_str());
						b_point = atoi(splitline[2].c_str());
						c_point = atoi(splitline[3].c_str());
					}
					//Make the points of triangle a,b,c.
					a = vert_list[a_point];
					b = vert_list[b_point];
					c = vert_list[c_point];

					glm::vec3 a_norm, b_norm, c_norm;
					if (OBJ_NORM) {
						a_norm = norm_list[a_point_norm];
						b_norm = norm_list[b_point_norm];
						c_norm = norm_list[c_point_norm];
					} else {
						a_norm.x = 0,a_norm.y = 0,a_norm.z = 0;
						b_norm.x = 0,b_norm.y = 0,b_norm.z = 0;
						c_norm.x = 0,c_norm.y = 0,c_norm.z = 0;
					}

					// Assuming that we can just use the matrix stack, but leave it as id_mat
					Transformation tri_trans(mat_stack);

					// Vertices of triangle.
					glm::vec3 vert_1 = tri_trans.world_point(vert_list[a_point]);
					glm::vec3 vert_2 = tri_trans.world_point(vert_list[b_point]);
					glm::vec3 vert_3 = tri_trans.world_point(vert_list[c_point]);

					// Normals of triangle.
					glm::vec3 norm_1 = tri_trans.world_normal(norm_list[a_point]);
					glm::vec3 norm_2 = tri_trans.world_normal(norm_list[b_point]);
					glm::vec3 norm_3 = tri_trans.world_normal(norm_list[c_point]);

					// Make pointer to this new triangle.
					// Hard coding BRDF
					ka.x = .1; ka.y = .3; ka.z = .9;
					kd.x = .3; kd.y = .3; kd.z = .9;
					ks.x = .3; kd.y = .3; kd.z = .5;
					kr.x = .2; kr.y = .2; kr.z = .2;
					float ior = 1.33;//if it is water
					ka.x = 0; ka.y = 0; ka.z = 0;
					kd.x = 0; kd.y = 0; kd.z = 0;
					ks.x = 0; kd.y = 0; kd.z = 0;
					kr.x = 0; kr.y = 0; kr.z = 0;
					sp = 30;
					RTriangle *t = new RTriangle(vert_1,vert_2,vert_3,ka,kd,ks,kr,ke,sp,norm_1,norm_2,norm_3);
					//Add triangle to scene.
					t->index_of_refraction = ior;
					t->transparency = true;//true;//not really working for triangles
					s.add_shape(t);

					// temp test
					if (flag) {
						float x = 0,y = 0, z = 0, r = 2;
						Sphere* sph = new Sphere(glm::vec3(x,y,z),r,ka,kd,ks,kr,ke,sp,tri_trans);
						//s.add_shape(sph);
						flag = false;
					}

					//resize vector. Obj files don't necessarily tell you how many faces there are going to be.
					// float n = glm::max(a_point,glm::max(b_point,c_point));
					// if(n+1>connected_triangles.size()){
					// 	connected_triangles.resize(n+1); //+1 because numbering for obj starts at 1.
					// }
				}
				else {
					std::cout << "Unknown command: " << splitline[0] << std::endl;
				}

			}
		}
	}
	// End Arg Parser

	// Set camera position if not specified (OBJ files)
	if (OBJ_ON) {
		// set_camera_and_perspective(c,max,min);

		// working camera basics, later set auto function
		float from_x = .25;
		float from_y = .25;
		float from_z = .65;
		float to_x = 0.25;
		float to_y = 0.25;
		float to_z = 0;
		float up_x = 0;
		float up_y = 1;
		float up_z = 0;
		float fov = 90;

		Camera cam(glm::vec3(from_x,from_y,from_z),glm::vec3(to_x,to_y,to_z),glm::vec3(up_x,up_y,up_z),fov);
		c = cam;


		//make checker board background
		RTriangle* t;
		glm::vec3 corner1(0,.5,0);
		glm::vec3 corner2(0,.25,0);
		glm::vec3 corner3(0,0,0);
		glm::vec3 corner4(.25,0,0);
		glm::vec3 corner5(.5,0,0);
		glm::vec3 corner6(.5,.25,0);
		glm::vec3 corner7(.5,.5,0);
		glm::vec3 corner8(.25,.5,0);
		glm::vec3 corner9(.25,.25,0);
		glm::vec3 ka(0,0,0); glm::vec3 d(0,0,0); glm::vec3 spec(0,0,0);
		glm::vec3 r(0,0,0);
		float sp = 0;

		glm::vec3 e(0,0,1);
		t = new RTriangle(corner1,corner2,corner9,ka,d,spec,r,e,sp);
		t->transparency = false;
		s.add_shape(t);
		e = glm::vec3(1,0,0);
		t = new RTriangle(corner2,corner3,corner9,ka,d,spec,r,e,sp);
		t->transparency = false;
		s.add_shape(t);
		e = glm::vec3(0,0,1);
		t->transparency = false;
		t = new RTriangle(corner3,corner4,corner9,ka,d,spec,r,e,sp);
		s.add_shape(t);
		e = glm::vec3(1,0,0);
		t->transparency = false;
		t = new RTriangle(corner4,corner5,corner9,ka,d,spec,r,e,sp);
		s.add_shape(t);
		e = glm::vec3(0,0,1);
		t->transparency = false;
		t = new RTriangle(corner5,corner6,corner9,ka,d,spec,r,e,sp);
		s.add_shape(t);
		e = glm::vec3(1,0,0);
		t->transparency = false;
		t = new RTriangle(corner6,corner7,corner9,ka,d,spec,r,e,sp);
		s.add_shape(t);
		e = glm::vec3(0,0,1);
		t->transparency = false;
		t = new RTriangle(corner7,corner8,corner9,ka,d,spec,r,e,sp);
		s.add_shape(t);
		e = glm::vec3(1,0,0);
		t->transparency = false;
		t = new RTriangle(corner8,corner1,corner9,ka,d,spec,r,e,sp);
		s.add_shape(t);
		t->transparency = false;

		//r = glm::vec3(1,1,1);
		//e = glm::vec3(0,0,0);
		//t = new RTriangle(glm::vec3(0,0,.1),glm::vec3(.5,0,.1),glm::vec3(.25,.4,.1),ka,d,spec,r,e,sp);
		//t->transparency = true;
		//t->index_of_refraction = 1.33;
		//s.add_shape(t);

		//t = new RTriangle(glm::vec3(0,0,.3),glm::vec3(.5,0,.3),glm::vec3(.25,.45,.3),ka,d,spec,r,e,sp);
		//t->transparency = true;
		//t->index_of_refraction = 1.33;
		//s.add_shape(t);


		//Transformation tri_trans(mat_stack);
		//Sphere* sph = new Sphere(glm::vec3(.25,.35,.2),.15,ka,d,spec,r,e,sp,tri_trans);
		//sph->index_of_refraction = 1.33;
		//sph->transparency = true;
		//s.add_shape(sph);
	}

	int BitsPerPixel = 24;
	Film canvas = Film(WIDTH, HEIGHT, BitsPerPixel, output_name);

	glm::vec3 UL,UR,LL,LR;

	c.cornerVectors(&UL,&UR,&LL,&LR,WIDTH,HEIGHT);
	s.set_params(c.position,UL,UR,LL,LR,WIDTH,HEIGHT,maxdepth);
	s.render(c,canvas);

	return 0;
}