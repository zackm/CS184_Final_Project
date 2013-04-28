#include "Scene.h"

#pragma once
#include <iostream>

#pragma once
#include <limits>

#pragma once
#include "BRDF.h"

#include "Sphere.h"

#include "Triangle.h"

#include "Shape.h"

#include "DirectionalLight.h"

using namespace std;

Scene::Scene(glm::vec3 eye,glm::vec3 UL_arg,glm::vec3 UR_arg, glm::vec3 LL_arg,
			 glm::vec3 LR_arg, int w,int h,int d) {

				 eye_position = eye;
				 UL = UL_arg;
				 UR = UR_arg;
				 LL = LL_arg;
				 LR = LR_arg;
				 width = w;
				 height = h;
				 maxdepth = d;
}

/*
For setting paramters after creating Scene object.
*/
void Scene::set_params(glm::vec3 eye,glm::vec3 UL_arg,glm::vec3 UR_arg, glm::vec3 LL_arg,
					   glm::vec3 LR_arg, int w,int h,int d) {

						   eye_position = eye;
						   UL = UL_arg;
						   UR = UR_arg;
						   LL = LL_arg;
						   LR = LR_arg;
						   width = w;
						   height = h;
						   maxdepth = d;
}

/*
Function iterates through each pixel and calls trace
to calculate the color of the pixel. We then commit 
the color to the Film kodak and generate the .png image.
*/
void Scene::render(Camera c, Film kodak) {
	float u,v;
	Ray ray;
	glm::vec3 color, pix_pos;

	float chunk = width / 100; // for progress indicator
	float counter = width / 100; // for progress indicator
	for (int i = 0; i < width; i++) {
		// for progress indicator
		if (i/1.0f > counter) {
			system("clear");
			cout << float(float(i) / width) * 100 << "% " << endl;
			counter += chunk;
		}

		for (int j = 0; j < height; j++) {

			//get u,v coordinates using the middle of the pixel
			u = float(float(i+.5)/width);
			v = float(float(j+.5)/height);
			pix_pos = u*(v*LL+(1-v)*UL)+(1-u)*(v*LR+(1-v)*UR);

			c.generateRay(pix_pos, &ray, eye_position);
			color = glm::vec3(0,0,0);

			trace(ray,&color);
			kodak.commit(width-i, height-j, color);
		}
	}
	system("clear");
	cout << "DONE" << endl;
	kodak.writeImage();
}

/*
Trace ray through scene. We use a while loop instead of a
purely recurisve call.
*/
void Scene::trace(Ray &r, glm::vec3 *color) {
	int i = 0;
	glm::vec3 reflection_coef(1.0f,1.0f,1.0f);
	float reflect_norm = 1.0;//to check if we need to break while loop early

	while (i<=maxdepth && reflect_norm>0.0f){
		float thit = std::numeric_limits<float>::infinity();
		LocalGeo local;
		Shape* best_shape;
		BRDF brdf;
		bool no_hit = true;
		glm::vec3 view_pos;
		for (std::list<Shape*>::iterator iter=shapes.begin(); iter != shapes.end(); ++iter) {
			Shape* s = *iter;

			float current_T;
			LocalGeo current_local;
			bool hit = (*s).intersect(r, &current_T, &current_local);
			if (current_T < thit && hit) {
				thit = current_T;
				local = current_local;
				no_hit = false;
				best_shape = s;
				view_pos = r.direction;
				brdf = s->get_brdf();
			}
		}

		if (no_hit) {
			//should break loop
			break;
		}

		if (i<1){
			//added only for initial rays
			*color += brdf.ka;
		}

		//added in all reflections
		*color += brdf.ke;

		Ray lray;
		glm::vec3 lcolor(0.0f,0.0f,0.0f);
		glm::vec3 add_color;

		for (std::list<Light*>::iterator iter=lights.begin(); iter != lights.end(); ++iter) {
			Light* l = *iter;

			(*l).generateLightRay(local,&lray,&lcolor);

			if (!intersect_checker(lray)) {
				add_color = shading(local, brdf, lray, lcolor, view_pos);
				*color += reflection_coef*add_color;
			}
		}

		//generate reflection ray and repeat
		//do this by resetting r.
		generateReflectionRay(local,&r);
		reflection_coef *= brdf.kr;
		reflect_norm = glm::dot(reflection_coef,reflection_coef);

		i++;
	}
}

/*
Updates ray arg so that we can reuse it inside of trace method.
*/
void Scene::generateReflectionRay(LocalGeo &local,Ray* ray){
	glm::vec3 normal = local.normal;
	glm::vec3 direction = -ray->direction;

	float norm_mag = glm::dot(normal,normal);
	if (norm_mag>0.0f){
		normal /= glm::sqrt(norm_mag);
	}

	float direction_mag = glm::dot(direction,direction);
	if (direction_mag>0.0f){
		direction /= glm::sqrt(direction_mag);
	}

	glm::vec3 reflection = (-direction)+2*glm::dot(direction,normal)*normal;
	float reflect_mag = glm::dot(reflection,reflection);
	if (reflect_mag>0.0f){
		reflection /= glm::sqrt(reflect_mag);
	}

	ray->direction = reflection;
	ray->direction = local.point;
	ray->t_min = .001;
	ray->t_max = std::numeric_limits<float>::infinity();
}

/*
Iterates through all shapes and just does basic intersect or not routine
*/
bool Scene::intersect_checker(Ray& r){
	for (std::list<Shape*>::iterator iter=shapes.begin(); iter != shapes.end(); ++iter) {
		Shape* s =  *iter;
		if ((*s).intersect(r)) {
			return true;
		}
	}
	return false;
}

/*
Using simple phong shading model.
*/
glm::vec3 Scene::shading(LocalGeo local, BRDF brdf, Ray lray, glm::vec3 lcolor,glm::vec3 view_pos){
	glm::vec3 normal = local.normal;
	glm::vec3 light_direction = lray.direction;

	float normal_mag = glm::sqrt(glm::dot(normal,normal));
	float light_mag = glm::sqrt(glm::dot(light_direction,light_direction));

	if (normal_mag>0){
		normal /= normal_mag;
	}
	if (light_mag>0){
		light_direction /= light_mag;
	}

	//Calculate the diffuse component
	float diffuse = glm::dot(normal,light_direction);

	glm::vec3 r_vec = (-light_direction)+(2.0f*diffuse)*normal;
	float r_norm = glm::dot(r_vec,r_vec);
	if (r_norm>0){
		r_vec /= glm::sqrt(r_norm);
	}

	diffuse = glm::max(diffuse,0.0f);

	//Calculate the specular component
	glm::vec3 view_vec = view_pos - local.point;
	float view_norm = glm::dot(view_vec,view_vec);
	if (view_norm > 0.0f) {
		view_vec /= glm::sqrt(view_norm);
	}
	float specular = glm::dot(r_vec,view_vec);
	specular = glm::max(specular,0.0f);
	specular = glm::pow(specular,brdf.shiny);

	glm::vec3 out_color;
	out_color[0] = (brdf.kd[0]*diffuse+brdf.ks[0]*specular)*lcolor[0];
	out_color[1] = (brdf.kd[1]*diffuse+brdf.ks[1]*specular)*lcolor[1];
	out_color[2] = (brdf.kd[2]*diffuse+brdf.ks[2]*specular)*lcolor[2];
	return out_color;
}

void Scene::add_shape(Shape* s) {
	shapes.push_front(s);
}

void Scene::add_light(Light* l) {
	lights.push_front(l);
}
