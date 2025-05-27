#pragma once
#include "../ThirdParty/pbrv4/pixelsensor.h"
#include "../ThirdParty/pbrv4/filters.h"
//holds information about the camera film: its dimensions, resolution, pixelsensor rgb color sensitivity, filter

struct pixel {
	glm::vec3 rgbsum = { 0.0f,0.0f,0.0f };
	float weightsum = 0;
};

class Film
{
public:
	std::vector<pixel> pixels;
	glm::ivec2 film_dim;
	glm::ivec2 image_res;
	pbrt::PixelSensor* pixel_sensor;
	pbrt::Filter* filter;
	//pdfs
};