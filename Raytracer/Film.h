#include "../FreeImage.h"
#include "glm/glm.hpp"

#pragma once
#include <string>

class Film{
public:
	FIBITMAP * bitmap;
	int BitsPerPixel;
	std::string filename;

	Film(int,int,int,std::string);
	void commit(int,int,glm::vec3);
	void writeImage();
};