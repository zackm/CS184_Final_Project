/* Film class that handles image output. */

#ifndef __FILM_H__
#define __FILM_H__
#include "Film.h"
#endif

#include "../FreeImage.h"
#include <iostream>
#include "glm/glm.hpp"

#pragma once
#include <string>

using namespace std;

Film::Film(int width, int height, int BitsPerPixel, std::string fname) {
	FreeImage_Initialise();
	bitmap = FreeImage_Allocate(width, height, BitsPerPixel);
	filename = fname;
}

void Film::commit(int x, int y, glm::vec3 color) {
	RGBQUAD free_image_color;

	//Cap color at max of 1.0f
	if (color[0] > 1.0f) {
		color[0] = 1.0f;
	}
	if (color[1] > 1.0f) {
		color[1] = 1.0f;
	}
	if (color[2] > 1.0f) {
		color[2] = 1.0f;
	}
	
	free_image_color.rgbRed = color[0] * 255.0f;
	free_image_color.rgbGreen = color[1] * 255.0f;
	free_image_color.rgbBlue = color[2] * 255.0f;
	FreeImage_SetPixelColor(bitmap, x, y, &free_image_color);

	return;
}

void Film::writeImage() {
	
	if (filename.empty()) {
		filename = "output.png";
	}

	if (FreeImage_Save(FIF_PNG, bitmap, filename.c_str(), 0)) 
		cout << "Image successfully saved!" << endl;

	FreeImage_DeInitialise();

	return;
}