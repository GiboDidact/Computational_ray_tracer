#pragma once
#include "../pch.h"
#include "Shapes.h"
#include "../ThirdParty/pbrv4/samplers.h"
/*

general:
every camera has to be able to go from renderspace to the film space to NDC and all those conversion
then it has to be able to generate ray samples at each pixel
has to also have its camera position and rotation transformations, and specific transformations

also shudder speed ill think on how to do that


idealized orthographic:
near/far plane, rays are just straight parralel lines from the pixels, not much just have to convert to NDC



idealized perspective:
near/far plane, you have a dot at the position and the near plane, the rays are just lines from the pinhole to the near plane, and you 
convert z values to NDC


pinhole:
radius of pinhole
dimension of the box
shuder?

calculating the light ray isn't too hard, you just pick a point on the pinhole and its just a straight line from that point to the point on the sensor


spherical:
ray directions are just the positions on the unit sphere in every direction
the only thing to really do is map the sphere to a plane, I can do spherical coordinates, or stereographic projection, other ones


thin lens physical model:
lens info: radius of curvature/focus point, lens diameter
apeture, size of the apeture maybe location or I can just put it practically in the lens
then the distance of the sensor, it can just be a variable

Finding the direction isn't too hard, one way is to just draw a line from sensor point to any point on lens and do snells law to refract
The other way without refracting uses the fact of the front focal plane, you just draw a line from the point on the sensor through the middle of the lens
then see where it intersects the focal plane, then thats the direction at the point of the lens that snells law would give

thick multilens physical model:



telescope/microscope model:

*/



/*
	1. create transformations between sample space to NDC to raster space, therefore I can map the sensor spots to pixels
	2. generate the samples based off location and direction
	3. given a location and direction calculate the outgoing ray origin and direction
	4. need the world space coordinates and perhaps transformations to camera space when I do the ray/object collision
	5. deal with extra camera stuff, focal lengths, apetures, shudder speeds, anything else
*/

//we can assume the little sensors on the film are 1 by 1 that way if its 100x100 each 1x1 is a sensor, you could give a difference size
//and really its just rescaling it

//I want to be able to control the sensor pixel resolution and size, that way I can control the data size maybe, I guess
//if you make sensor pixels bigger you'd need more sample points so it wouldn't really save time and it would all be blurred together

//camera sensors dont pick up a spectral distribution data, it just picks up 3 scalar RGB color space, though that space can be different than the output color
//space of the camera

//Screen space is the film in world space, z usually between 0 and 1
//NDC (0,0) is topleft normalized values
//Raster (0,0) is topleft but scaled
class CameraBase
{
public:
	CameraBase(float sensor_width, float sensor_height, glm::vec3 wrld_position, glm::vec3 look_d, glm::vec3 r_dir,
		glm::vec3 world_updir, glm::vec2 image_res, float _lensRadius = 0, float _focalDistance = 0) : sensor_dimensions(sensor_width, sensor_height), world_pos(wrld_position),
		look_direction(look_d), right_direction(r_dir), worldup_direction(world_updir), image_resolution(image_res), M_CameratoScreen(glm::mat4(1.0f)),
		M_RastertoCamera(glm::mat4(1.0f)), lensRadius(_lensRadius), focalDistance(_focalDistance)
	{
		up_direction = glm::cross(look_direction, right_direction);

		//screen-->raster
		glm::mat4 M_ScreentoNDC = glm::scale(glm::mat4(1.0f), glm::vec3(1.0f / sensor_dimensions.x, 1.0f / sensor_dimensions.y, 1))
			* glm::translate(glm::mat4(1.0f), glm::vec3(sensor_dimensions.x / 2.0f, sensor_dimensions.y / 2.0f, 0));
		glm::mat4 M_NDCtoRaster = glm::scale(glm::mat4(1.0f), glm::vec3(image_resolution.x, -image_resolution.y, 1))
			* glm::translate(glm::mat4(1.0f), glm::vec3(0, -1, 0));
		
		M_ScreentoRaster = M_NDCtoRaster * M_ScreentoNDC;
		M_RastertoScreen = glm::inverse(M_ScreentoRaster);
	
		calculateWorldCameraMatrices();

		if (false)
		{
			std::cout << "updirection_y: "<< up_direction.y << std::endl;
			
			//screen to raster
			glm::vec4 _point1(250.0f, 250.0f, 0, 1); //map to (0,600)
			glm::vec4 point1 = (M_ScreentoRaster) * _point1;
			std::cout << point1.x << " " << point1.y << " " << point1.z << " " << point1.z << std::endl;

			//raster to screen
			glm::vec4 _point2(0, 0, 0, 1); //map to (0,600)
			glm::vec4 point2 = (M_RastertoScreen)*_point2;
			std::cout << point2.x << " " << point2.y << " " << point2.z << " " << point2.z << std::endl;
		}
	}
	~CameraBase() = default;

	void setWorldPos(glm::vec3 pos)
	{
		world_pos = pos;
		calculateWorldCameraMatrices();
	}

	void setyawpitch(float yaw, float pitch)
	{
		look_direction.x = sin(glm::radians(pitch)) * cos(glm::radians(yaw));
		look_direction.y = cos(glm::radians(pitch));
		look_direction.z = sin(glm::radians(pitch)) * sin(glm::radians(yaw));
		look_direction = glm::normalize(look_direction);
		calculateWorldCameraMatrices();
	}

	void calculateWorldCameraMatrices()
	{
		glm::vec3 dir = glm::normalize(look_direction);
		right_direction = glm::normalize(glm::cross(worldup_direction, dir));
		up_direction = glm::cross(dir, right_direction);

		M_CameratoWorld = glm::mat4(glm::vec4(right_direction.x, right_direction.y, right_direction.z, 0),
			glm::vec4(up_direction.x, up_direction.y, up_direction.z, 0),
			glm::vec4(dir.x, dir.y, dir.z, 0),
			glm::vec4(world_pos.x, world_pos.y, world_pos.z, 1));

		M_WorldtoCamera = glm::inverse(M_CameratoWorld);


		if (false)
		{
			std::cout << "world_pos: " << world_pos.x << " " << world_pos.y << " " << world_pos.z << std::endl;
			std::cout << "dir: " << dir.x << " " << dir.y << " " << dir.z << std::endl;
			std::cout << "right: " << right_direction.x << " " << right_direction.y << " " << right_direction.z << std::endl;
			std::cout << "newup: " << up_direction.x << " " << up_direction.y << " " << up_direction.z << std::endl;

			glm::vec4 _pos1(5, 5, 5, 1);
			glm::vec4 pos1 = M_CameratoWorld * _pos1;
			std::cout << pos1.x << " " << pos1.y << " " << pos1.z << " " << pos1.w << std::endl;
			char a;
			std::cin >> a;
		}
	}
	
	glm::mat4 GetWorldtoCameraMatrix()
	{
		return M_WorldtoCamera;
	}

	glm::mat4 GetCameratoWorldMatrix()
	{
		return M_CameratoWorld;
	}

	glm::mat4 GetWorldtoCameraWorldMatrix()
	{
		//converts the position of camera, but keeps world orientation
		return glm::translate(glm::mat4(1.0f), -world_pos);
	}

	glm::vec3 getlookdirection() { return look_direction; }
	glm::vec3 getrightdirection() { return right_direction; }
	glm::vec3 getupdirection() { return up_direction; }

	virtual Ray generateRay(glm::vec2 pixel, pbrt::Sampler* sampler) = 0;

	void SetlensRadius(float r) { lensRadius = r; }
	void SetfocalDistance(float d) { focalDistance = d; }
	//camera transform
	//film
	//shutter
	//medium
	//differentialXY

protected:
	float lensRadius;
	float focalDistance;
	//normalize directions
	glm::vec3 world_pos; //(x,y,z)
	glm::vec3 look_direction;
	glm::vec3 right_direction;
	glm::vec3 up_direction;
	glm::vec3 worldup_direction; //helps us recalculate right_direction

	glm::vec2 sensor_dimensions; //(width,height)
	glm::vec2 image_resolution;

	//transforms
	glm::mat4 M_ScreentoRaster;
	glm::mat4 M_RastertoScreen;

	glm::mat4 M_CameratoWorld;
	glm::mat4 M_WorldtoCamera;

	glm::mat4 M_CameratoScreen;
	glm::mat4 M_RastertoCamera;
};

class OrthographicCamera : public CameraBase
{
public:
	OrthographicCamera(float near, float far, float sensor_width, float sensor_height, glm::vec3 wrld_position, glm::vec3 look_d, glm::vec3 r_dir,
		glm::vec3 world_updir, glm::vec2 image_res) : CameraBase(sensor_width, sensor_height, wrld_position, look_d, r_dir, world_updir, image_res)
	{
		N = near;
		F = far;
		
		//orthographic projection, x,y is same, z is just scaled from near to far (0,1)
		M_CameratoScreen = glm::scale(glm::mat4(1.0f), glm::vec3(1, 1, 1.0 / (F - N)))
			               * glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -N));
		
		M_RastertoCamera = glm::inverse(M_CameratoScreen) * M_RastertoScreen;
	}
	~OrthographicCamera() = default;

	Ray generateRay(glm::vec2 pixel, pbrt::Sampler* sampler) override
	{
		//has to have w=1
		glm::vec4 Cam_pos = M_RastertoCamera * glm::vec4(pixel.x, pixel.y, 0, 1);
		if (Cam_pos.w != 1)
			std::cout << "orthoray w not 1\n";
		Ray ray(Cam_pos, glm::vec3(0, 0, 1));

		//convert to world space
		ray.Transform(M_CameratoWorld);

		return ray;
	}
private:
	float N, F;	
};
//Near= 1e-2f, Far=1000f

class PerspectiveCamera : public CameraBase
{
public:
	//you still have a sensor_width/height thats the near plane, the fov will then just truncate some of
	//its just scaling like 2 dof with 3 variables: sensor_width/height, fov, near plane
	PerspectiveCamera(float near, float far, float sensor_width, float sensor_height, float _fov, glm::vec3 wrld_position, glm::vec3 look_d, 
		glm::vec3 r_dir, glm::vec3 world_updir, glm::vec2 image_res, float _lensRadius = 0, float _focalDistance = 0) :
		CameraBase(2 * near * tan(glm::radians(_fov) / 2.0f), 2 * near * tan(glm::radians(_fov) / 2.0f) * (image_res.x / (float)image_res.y), wrld_position, look_d, r_dir, world_updir, image_res, _lensRadius, _focalDistance)
	{
		//std::cout << 2*near*std::tan(glm::radians(_fov) / 2.0f) << std::endl;

		N = near;
		F = far;
		fov = _fov;

		calculuateMatrix();
	}
	~PerspectiveCamera() = default;

	void ChangeFOV(float _fov)
	{
		fov = _fov;
		calculuateMatrix();
	}

	Ray generateRay(glm::vec2 pixel, pbrt::Sampler* sampler) override
	{
		//divide by w
		glm::vec4 near_pos_w = M_RastertoCamera * glm::vec4(pixel.x, pixel.y, 0, 1);
		glm::vec3 near_pos = glm::vec3(near_pos_w.x / near_pos_w.w, near_pos_w.y / near_pos_w.w, near_pos_w.z / near_pos_w.w);
		
		Ray ray(glm::vec3(0, 0, 0), glm::normalize(near_pos));

		if (lensRadius > 0 && sampler)
		{
			glm::vec2 lens_pos = glm::vec2(lensRadius, lensRadius) * SampleUniformDiskConcentric(sampler->Get2D());

			//really focaldistance is the focal plane, not the near plane, so you want to find the point on the focal plane and subtract that from lenspos
			float ft = focalDistance / ray.d.z;
			glm::vec3 pfocus = ray.o + ray.d * ft;

			ray.o = glm::vec3(lens_pos.x, lens_pos.y, 0);
			ray.d = glm::normalize(pfocus - ray.o);
		}

		//convert to world space
		ray.Transform(M_CameratoWorld);

		return ray;
	}

private:
	float N, F;
	float fov;

	void calculuateMatrix()
	{
		float invTanAng = 1.0f / std::tan(glm::radians(fov) / 2.0f);
		glm::mat4 perspective = glm::mat4(glm::vec4(1, 0, 0, 0), glm::vec4(0, 1, 0, 0), glm::vec4(0, 0, F / (F - N), 1), glm::vec4(0, 0, -(F * N / (F - N)), 0));
		M_CameratoScreen = glm::scale(perspective, glm::vec3(invTanAng, invTanAng, 1));

		M_RastertoCamera = glm::inverse(M_CameratoScreen) * M_RastertoScreen;
	}
};

class PinholeCamera : public CameraBase
{
public:
	PinholeCamera(float radius, glm::vec3 _box_dimensions, glm::vec3 wrld_position, glm::vec3 look_d, glm::vec3 r_dir,
		glm::vec3 world_updir, glm::vec2 image_res) : CameraBase(_box_dimensions.x, _box_dimensions.y, wrld_position, look_d, r_dir, world_updir, image_res)
	{
		hole_radius = radius;
		box_dimensions = _box_dimensions;

		//don't need these
		M_CameratoScreen = glm::mat4(1.0f);
		M_RastertoCamera = glm::mat4(1.0f);
	}
	~PinholeCamera() = default;

	Ray generateRay(glm::vec2 pixel, pbrt::Sampler* sampler) override
	{
		glm::vec3 sensor_pos = M_RastertoScreen * glm::vec4(pixel.x, pixel.y, 0, 1);
		//the pinhole position could be anything within the hole_radius
		glm::vec3 pinhole_pos = glm::vec3(0.0f * hole_radius * cos(glm::radians(0.0f)), 0.0f * hole_radius * sin(glm::radians(0.0f)),
			box_dimensions.z);
		Ray ray(sensor_pos, glm::normalize(pinhole_pos - sensor_pos));

		//convert to world space
		ray.Transform(M_CameratoWorld);

		return ray;
	}

	Ray generateRay(glm::vec2 pixel, float lens_angle, float len_percent_r)
	{
		glm::vec3 sensor_pos = M_RastertoScreen*glm::vec4(pixel.x, pixel.y, 0, 1);
		//the pinhole position could be anything within the hole_radius
		glm::vec3 pinhole_pos = glm::vec3(len_percent_r*hole_radius*cos(glm::radians(lens_angle)), len_percent_r*hole_radius* sin(glm::radians(lens_angle)), 
			                              box_dimensions.z);
		Ray ray(sensor_pos, glm::normalize(pinhole_pos - sensor_pos));

		//convert to world space
		ray.Transform(M_CameratoWorld);

		return ray;
	}

private:
	float hole_radius;
	glm::vec3 box_dimensions; //(w,h,l)
};

//1 positive lens, with an apeture
class ThinlensCamera : public CameraBase
{
public:
	ThinlensCamera(float Radius_curvature, float lens_d, float _apeture, float s_depth, float sensor_width, float sensor_height, glm::vec3 wrld_position, glm::vec3 look_d,
		glm::vec3 r_dir, glm::vec3 world_updir, glm::vec2 image_res) : CameraBase(sensor_width, sensor_height, wrld_position, look_d, r_dir, world_updir, image_res)
	{
		R = Radius_curvature;
		F = R / 2.0f; //focal point can be calculated from radius of curvature 
		lens_diameter = lens_d;
		
		apeture = _apeture;
		
		sensor_depth = s_depth;
	}
	~ThinlensCamera() = default;

	Ray generateRay(glm::vec2 pixel, float lens_angle, float len_percent_r)
	{
		//The other way without refracting uses the fact of the front focal plane, you just draw a line from the point on the sensor through the middle of the lens
		//then see where it intersects the focal plane, then thats the direction at the point of the lens that snells law would give
		
		glm::vec3 sensor_pos = M_RastertoScreen * glm::vec4(pixel.x, pixel.y, 0, 1);
		//with apeture the lense diameter is just smaller
		float apeture_diameter = lens_diameter - apeture;
		glm::vec3 lens_pos = glm::vec3(len_percent_r * (apeture_diameter / 2.0f) * cos(glm::radians(lens_angle)),
			len_percent_r * (apeture_diameter / 2.0f) * sin(glm::radians(lens_angle)), sensor_depth);
		glm::vec3 lens_center = glm::vec3(0, 0, sensor_depth);

		glm::vec3 tmpray = glm::normalize(lens_center - sensor_pos);
		float t = F / tmpray.z;
		glm::vec3 focal_plane_pos = lens_center + tmpray * t;
		Ray ray(lens_pos, glm::normalize(focal_plane_pos - lens_pos));
		
		//convert to world space
		ray.Transform(M_CameratoWorld);

		return ray;
	}

private:
	float R;
	float F; //on z axis
	float lens_diameter;

	float apeture;
	
	float sensor_depth;
};

class SphericalCamera : public CameraBase
{

};

class ThicklensCamera : public CameraBase
{

};

class TelescopeCamera : public CameraBase
{

};