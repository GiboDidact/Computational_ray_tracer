#pragma once
#include "../pch.h"
#include "../ThirdParty/pbrv4/helpers.h"
#include "AssetManager.h"

/*

each shape needs a ray/intersection function that gives the collision point
then at each point: (u,v), partial derivatives, the normal, the tangent plane

For now I just need: position hit, (u,v) values, and the normal 
normal is just normalize(cross(t_u,t_v)), or if its an implicit equation its the gradient of f(x,y,z) = k 


sphere
cylinder
disks
triangle meshes
bilinear patches
curves
subdivision mesh

*/

//objects you pass in its world transformation (rotating,scaling,translating), and going backwards from world to object space is inverse of that

//make sure transformations and normals and rays are preserved, keep that in mind

//their object space z is up, even though z is forward in world so whatever I guess you can always rotate it

//*might have to deal with self-intersection, ray tracing gems has a chapter on it

struct Ray
{
	Ray(glm::vec3 origin, glm::vec3 direction) : o(origin), d(direction) { }

	void Transform(glm::mat4 M)
	{
		o = M * glm::vec4(o.x, o.y, o.z, 1);
		d = glm::normalize(M * glm::vec4(d.x, d.y, d.z, 0));
	}

	Ray operator*(const glm::mat4& mat)
	{
		Ray new_ray(mat*glm::vec4(o.x,o.y,o.z,1), mat * glm::vec4(d.x, d.y, d.z, 1));
	}

	glm::vec3 o;
	glm::vec3 d; //normalized
};


class Bounds3
{
public:
	Bounds3() : pmin(0, 0, 0), pmax(0, 0, 0) {}
	Bounds3(glm::vec3 min, glm::vec3 max) : pmin(min), pmax(max) {}
	~Bounds3() = default;

	void Transform(glm::mat4 M)
	{
		//move back to origin then transform then move back for scaling, dont think we need to do it here
		glm::vec3 half_d = glm::vec3(pmax.x - pmin.x, pmax.y - pmin.y, pmax.z - pmin.z) / 2.0f;
		glm::vec3 C = pmin + half_d;
		glm::mat4 T = glm::translate(glm::mat4(1.0f), -C);
		glm::mat4 T_I = glm::translate(glm::mat4(1.0f), C);
		//M = T_I * M * T;

		//scaling and translating you can just check min and max
		//its an axis-aligned box so we need to rotate all points and find min/max
		std::vector<glm::vec3> transformed_points = { M * glm::vec4(pmin.x, pmin.y, pmin.z, 1),
										  M * glm::vec4(pmin.x, pmax.y, pmin.z, 1),
										  M * glm::vec4(pmin.x, pmax.y, pmax.z, 1),
										  M * glm::vec4(pmin.x, pmin.y, pmax.z, 1),
										  M * glm::vec4(pmax.x, pmax.y, pmax.z, 1),
										  M * glm::vec4(pmax.x, pmin.y, pmax.z, 1),
										  M * glm::vec4(pmax.x, pmin.y, pmin.z, 1),
										  M * glm::vec4(pmax.x, pmax.y, pmin.z, 1) };
		
		glm::vec2 x_minmax(std::numeric_limits<float>::max(), std::numeric_limits<float>::min());
		glm::vec2 y_minmax(std::numeric_limits<float>::max(), std::numeric_limits<float>::min());
		glm::vec2 z_minmax(std::numeric_limits<float>::max(), std::numeric_limits<float>::min());
		for (int i = 0; i < transformed_points.size(); i++)
		{
			x_minmax.x = std::min(x_minmax.x, transformed_points[i].x);
			x_minmax.y = std::max(x_minmax.y, transformed_points[i].x);

			y_minmax.x = std::min(y_minmax.x, transformed_points[i].y);
			y_minmax.y = std::max(y_minmax.y, transformed_points[i].y);

			z_minmax.x = std::min(z_minmax.x, transformed_points[i].z);
			z_minmax.y = std::max(z_minmax.y, transformed_points[i].z);
		}


		pmin = glm::vec3(x_minmax.x, y_minmax.x, z_minmax.x);
		pmax = glm::vec3(x_minmax.y, y_minmax.y, z_minmax.y);
	}

	bool IntersectP(const Ray& ray, float tMax = std::numeric_limits<float>::max(), float* hitt0 = nullptr, float* hitt1 = nullptr) const
	{
		//get 6 intersection points, 2 for each dimension. then see if each intersection is inside the thing if it is you can early true
		//the equation for the ray to hit some plane on a dimension is just t = x/y/z- ray.x/y/z / direction.x/y/z
		//you get an interval of t's then any new intersection has to be within that if its outside well we know it cant be in the box

		float min_t = 0, max_t = tMax;
		for (int i = 0; i < 3; i++)
		{
			float invRayDir = 1 / ray.d[i];
			float tNear = (pmin[i] - ray.o[i]) * invRayDir;
			float tFar = (pmax[i] - ray.o[i]) * invRayDir;

			if (tNear > tFar) std::swap(tNear, tFar);
			tFar *= 1 + 2 * pbrt::gamma(3);

			min_t = tNear > min_t ? tNear : min_t;
			max_t = tFar < max_t ? tFar : max_t;
			if (min_t > max_t) return false;
		}
		if (hitt0) *hitt0 = min_t;
		if (hitt1) *hitt1 = max_t;

		return true;
	}

	glm::vec3 pmin, pmax;
};

inline Bounds3 TransformBounds(Bounds3 b, const glm::mat4& M)
{
	Bounds3 bounds = b;
	bounds.Transform(M);

	return bounds;
}
inline Ray TransformRay(Ray r, const glm::mat4& M)
{
	Ray ray = r;
	ray.Transform(M);

	return ray;
}

struct LocalSurfaceInfo
{
	//tranforms usually from object space to world space
	void Transform(glm::mat4 M)
	{
		//normals
		glm::mat3 normal_transform = glm::transpose(glm::inverse(M));
		n = glm::normalize(normal_transform * n);
		wo = glm::normalize(normal_transform * wo);

		//positions
		hitp = M * glm::vec4(hitp, 1.0f);
		
		//vectors
		du = M * glm::vec4(du, 0);
		dv = M * glm::vec4(dv, 0);
	}

	float tHit;
	glm::vec3 hitp;
	float u, v;
	glm::vec3 du, dv;
	glm::vec3 n;
	glm::vec3 wo;
	//material
	// partial normal derivatives
};

class Shape
{
public:
	Shape(const std::string& _name, glm::mat4 rigidtransform) : name(_name)
	{
		//world space is z forward, object space is z up, so z and y have to be flipped
		glm::mat4 permuation_y_z(1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1);
		//ObjectToRender = permuation_y_z * rigidtransform;
		ObjectToRender = rigidtransform * permuation_y_z;
		RenderToObject = glm::inverse(ObjectToRender);
	}
	~Shape() = default;

	virtual void SetRigidTransform(glm::mat4 rigidtransform)
	{
		glm::mat4 permuation_y_z(1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1);
		ObjectToRender = rigidtransform * permuation_y_z;
		RenderToObject = glm::inverse(ObjectToRender);
	}

	std::string GetName() const { return name; }
	glm::mat4 GetRenderToObjectMatrix() const { return RenderToObject; }
	glm::mat4 GetObjectToRenderMatrix() const { return ObjectToRender; }
	virtual Bounds3 Bounds() const = 0;
	virtual std::optional<LocalSurfaceInfo> Intersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const = 0;
	virtual bool IntersectP(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const = 0;
	virtual float Area() const = 0;
	
	//sample
	//pdf
	//normal bounds

protected:
	std::string name;
	glm::mat4 RenderToObject, ObjectToRender;
};

class Sphere : public Shape
{
public:
	struct SphereIntersect
	{
		float t;
		glm::vec3 hitp;
		glm::vec3 ray_d;
		float phi;
	};
public:
	Sphere(const std::string& _name, glm::mat4 rigidtransform, float radius, float _zmin, float _zmax, float _phimax) :
		Shape(_name, rigidtransform) 
	{ 
		r = radius;
		zmin = glm::clamp(_zmin, -r, r);
		zmax = glm::clamp(_zmax, -r, r);
		thetamin = std::acos(glm::clamp(zmin / r, -1.f, 1.f));
		thetamax = std::acos(glm::clamp(zmax / r, -1.f, 1.f));

		//std::cout << thetamin << " " << thetamax << std::endl;
		phimax = glm::radians(glm::clamp(_phimax, 0.0f, 360.f));
	}
	~Sphere() = default;

	float Area() const 
	{ 
		return phimax * r * (zmax - zmin);
	}
	
	Bounds3 Bounds() const 
	{
		return TransformBounds(Bounds3(glm::vec3(-r, -r, zmin), glm::vec3(r, r, zmax)), ObjectToRender);
	}

	std::optional<LocalSurfaceInfo> Intersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		std::optional<SphereIntersect> isect = BasicIntersect(ray, tMax);
		if (!isect)
			return {};

		glm::vec3 p = isect.value().hitp;

		//normal in world space, position in world space
		LocalSurfaceInfo info;
		info.tHit = isect.value().t;
		info.hitp = p;
		glm::vec2 uv = PostoUV(p.x, p.y, p.z);
		info.u = glm::clamp(uv.x,0.0f,1.0f);
		info.v = glm::clamp(uv.y, 0.0f, 1.0f);
		info.du = calculate_du(p.x, p.y, p.z);
		info.dv = calculate_dv(p.x, p.y, p.z);
		info.n = calculate_normal(p.x, p.y, p.z);
		if (glm::dot(info.n, isect.value().ray_d) > 0) //the normal and the ray direction should be on opposite sides
			info.n = -info.n;
		info.wo = glm::vec3(0, 0, 1);

		//transform to world space
		info.Transform(ObjectToRender);

		return LocalSurfaceInfo{ info };
	}

	bool IntersectP(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		return BasicIntersect(ray, tMax).has_value();
	}

	std::optional<SphereIntersect> BasicIntersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		//transform ray from rendertoobject space
		glm::vec3 o = RenderToObject * glm::vec4(ray.o.x, ray.o.y, ray.o.z, 1);
		glm::vec3 d = RenderToObject * glm::vec4(ray.d.x, ray.d.y, ray.d.z, 0);
		//std::cout << "o: " << o.x << " " << o.y << " " << o.z << std::endl;
		//std::cout << "d: " << d.x << " " << d.y << " " << d.z << std::endl;

		//solve the quadratic and get t0, t1 (from pbrv4, a more robust way to solve the quadratic)
		float a = d.x * d.x + d.y * d.y + d.z*d.z;
		float b = 2 * (d.x * o.x + d.y * o.y + d.z*o.z);
		float c = o.x * o.x + o.y * o.y + o.z*o.z - r * r;

		glm::vec3 v = o - b / (2 * a) * d;
		float length = glm::length(v);
		float discrim = 4 * a * (r + length) * (r - length);
		//std::cout << discrim << std::endl;
		if (discrim < 0)
			return {};

		float rootDiscrim = std::sqrt(discrim);
		float q;
		if (b < 0)
			q = -.5f * (b - rootDiscrim);
		else
			q = -.5f * (b + rootDiscrim);
		float t0 = q / a;
		float t1 = c / q;

		//choose the nearest hit and get intersection point
		if (t0 > t1)
		{
			float tmp = t1;
			t1 = t0;
			t0 = tmp;
		}
		if (t0 > tMax || t1 <= 0)
			return {};
		float tShapeHit = t0;
		//std::cout << t0 << " " << t1 << std::endl;
		if (tShapeHit <= 0)
		{
			tShapeHit = t1;
			if (tShapeHit > tMax)
				return {};
		}

		glm::vec3 hitp = o + tShapeHit * d;
		//std::cout << hitp.x << " " << hitp.y << std::endl;
		//refine sphere intersection point
		hitp *= r / glm::distance(hitp, glm::vec3(0, 0, 0));

		if (hitp.x == 0 && hitp.y == 0)
			hitp.x = 1e-5 * r;
		float phi = std::atan2(hitp.y, hitp.x);
		if (phi < 0) phi += 2 * std::numbers::pi;

		//check clipping
		if (hitp.z < zmin || hitp.z > zmax || phi > phimax)
		{
			//got clipped, see if other t1 will work
			if (tShapeHit == t1)
				return {};
			if (t1 > tMax)
				return {};
			tShapeHit = t1;

			hitp = o + tShapeHit * d;
			hitp *= r / glm::distance(hitp, glm::vec3(0, 0, 0));

			if (hitp.x == 0 && hitp.y == 0)
				hitp.x = 1e-5 * r;
			phi = std::atan2(hitp.y, hitp.x);
			if (phi < 0) phi += 2 * std::numbers::pi;

			if (hitp.z < zmin || hitp.z > zmax || phi > phimax)
				return {};
		}
		
		return SphereIntersect{tShapeHit, hitp, glm::normalize(d), phi};
	}

	glm::vec3 UVtoPos(float u, float v) const
	{
		float phi = u * phimax;
		float theta = thetamin + v * (thetamax - thetamin);
		glm::vec3 pos(r *sin(theta)*cos(phi), r * sin(theta)*sin(phi), r*cos(theta));

		return pos;
	}

	//assuming x,y,z is on the shape
	glm::vec2 PostoUV(float x, float y, float z) const
	{
		glm::vec2 angles = PostoAngles(x, y, z);
		float theta = angles.x;
		float phi = angles.y;
		
		float u = phi / phimax;
		float v = (theta - thetamin) / (thetamax - thetamin);

		return glm::vec2(u, v);
	}
	glm::vec2 PostoAngles(float x, float y, float z) const
	{
		float theta = acos(glm::clamp(z / r, -1.f, 1.f));
		float phi = std::atan2(y, x);
		if (phi < 0) phi += 2 * std::numbers::pi;

		return glm::vec2(theta, phi);
	}

	//glm::vec3 calculate_du(float u, float v) const
	//{
		//return glm::normalize(glm::vec3(-r * max_phi * sin(u * max_phi), r * max_phi * cos(u * max_phi), 0));
	//}
	glm::vec3 calculate_du(float x, float y, float z) const
	{
		return glm::normalize(glm::vec3(-phimax * y, phimax * x, 0));
	}
	//glm::vec3 calculate_dv(float u, float v) const
	//{
		//return glm::normalize(glm::vec3(0, 0, max_z - min_z));
	//}
	glm::vec3 calculate_dv(float x, float y, float z) const
	{
		glm::vec2 angles = PostoAngles(x, y, z);
		float theta = angles.x;
		float phi = angles.y;

		return glm::normalize((thetamax-thetamin)*glm::vec3(z*cos(phi), z*sin(phi), -r*sin(theta)));
	}

	glm::vec3 calculate_normal(float u, float v) const
	{
		glm::vec3 p = UVtoPos(u, v);
		return glm::normalize(glm::cross(calculate_du(p.x, p.y, p.z), calculate_dv(p.x, p.y, p.z)));
	}

	glm::vec3 calculate_normal(float x, float y, float z) const
	{
		//gradient f(x,y,z) = 0, (df/dx,df/dy,dy/dz)
		glm::vec3 N = glm::vec3(2 * x, 2 * y, 2 * z);
		return glm::normalize(N);
	}

private:
	//implicit function is x^2+y^2+z^2-r^2 = 0
	//parameterization is x=rsin(theta)cos(phi), y=rsin(theta)sin(phi), z=rcos(theta),  theta(0,pi) phi(0,2pi)
	//phi = u*phi_max
	//theta = thetamin + v*(thetamax-thetamin)

	float r;
	float zmin, zmax;
	float thetamin, thetamax, phimax;
};

class Cylinder : public Shape
{
public:
	struct CylinderIntersect
	{
		float t;
		glm::vec3 hitp;
		glm::vec3 ray_d;
		float phi;
	};
public:
	Cylinder(const std::string& _name, glm::mat4 rigidtransform, float radius, float z_min, float z_max, float phi_max)
		: Shape(_name, rigidtransform) 
	{ 
		r = radius;
		min_z = z_min;
		max_z = z_max;
		max_phi = glm::radians(phi_max);
	}
	~Cylinder() = default;

	float Area() const 
	{ 
		return  (max_z - min_z) * r * max_phi;
	}
	Bounds3 Bounds() const 
	{ 
		return TransformBounds(Bounds3(glm::vec3(-r, -r, min_z), glm::vec3(r, r, max_z)),ObjectToRender);
	}
	
	std::optional<LocalSurfaceInfo> Intersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		std::optional<CylinderIntersect> isect = BasicIntersect(ray, tMax);
		if (!isect)
			return {};

		glm::vec3 p = isect.value().hitp;

		LocalSurfaceInfo info;
		info.tHit = isect.value().t;
		info.hitp = p;
		glm::vec2 uv = PostoUV(p.x, p.y, p.z);
		info.u = glm::clamp(uv.x, 0.0f, 1.0f);
		info.v = glm::clamp(uv.y, 0.0f, 1.0f);
		info.du = calculate_du(p.x, p.y, p.z);
		info.dv = calculate_dv(p.x, p.y, p.z);
		info.n = calculate_normal(p.x, p.y, p.z);
		if (glm::dot(info.n, isect.value().ray_d) > 0) //the normal and the ray direction should be on opposite sides
		{
			info.n = -info.n;
		}

		info.wo = glm::vec3(0, 0, 1);

		//transform to world space
		info.Transform(ObjectToRender);

		return LocalSurfaceInfo{ info };
	}

	bool IntersectP(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		return BasicIntersect(ray, tMax).has_value();
	}

	std::optional<CylinderIntersect> BasicIntersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		//transform ray from rendertoobject space
		glm::vec3 o = RenderToObject * glm::vec4(ray.o.x, ray.o.y, ray.o.z,1);
		glm::vec3 d = RenderToObject * glm::vec4(ray.d.x, ray.d.y, ray.d.z, 0);

		//solve the quadratic and get t0, t1 (from pbrv4, a more robust way to solve the quadratic)
		float a = d.x * d.x + d.y * d.y;
		float b = 2 * (d.x * o.x + d.y * o.y);
		float c = o.x * o.x + o.y * o.y - r * r;

		float f = b / (2 * a);
		float vx = o.x - f * d.x, vy = o.y - f * d.y;
		float length = std::sqrtf(vx * vx + vy * vy);
		float discrim = 4 * a * (r + length) * (r - length);
		if (discrim < 0)
			return {};

		float rootDiscrim = std::sqrt(discrim);
		float q;
		if (b < 0)
			q = -.5f * (b - rootDiscrim);
		else
			q = -.5f * (b + rootDiscrim);
		float t0 = q / a;
		float t1 = c / q;

		//choose the nearest hit
		if (t0 > t1)
		{
			float tmp = t1;
			t1 = t0;
			t0 = tmp;
		}
		if (t0 > tMax || t1 <= 0)
			return {};
		float tShapeHit = t0;
		if (tShapeHit <= 0)
		{
			tShapeHit = t1;
			if (tShapeHit > tMax)
				return {};
		}

		//check to make sure its in the z range and the angle range
		glm::vec3 hit_p = o + tShapeHit * d;
		float phi = std::atan2(hit_p.y, hit_p.x);
		if (phi < 0) phi += 2 * std::numbers::pi;

		if (hit_p.z < min_z || hit_p.z > max_z || phi > max_phi)
		{
			//one of the t's failed, try the other if we can
			if (tShapeHit == t1)
				return {};
			tShapeHit = t1;
			if (t1 > tMax)
				return {};
			hit_p = o + tShapeHit * d;
			phi = std::atan2(hit_p.y, hit_p.x);
			if (phi < 0) phi += 2 * std::numbers::pi;
			if (hit_p.z < min_z || hit_p.z > max_z || phi > max_phi)
				return {};
		}

		return CylinderIntersect{tShapeHit, hit_p, d, phi };
	}

	glm::vec3 UVtoPos(float u, float v) const
	{
		float phi = u * max_phi;
		glm::vec3 pos(r * cos(phi), r * sin(phi), min_z + v * (max_z - min_z));
		
		return pos;
	}

	//assuming x,y,z is on the shape
	glm::vec2 PostoUV(float x, float y, float z) const
	{
		float phi = std::atan2(y, x);
		if (phi < 0) phi += 2 * std::numbers::pi;
		float u = phi / max_phi;
		float v = (z - min_z) / (max_z - min_z);
		
		return glm::vec2(u, v);
	}

	glm::vec3 calculate_du(float u, float v) const
	{
		return glm::normalize(glm::vec3(-r * max_phi * sin(u * max_phi), r * max_phi * cos(u * max_phi), 0));
	}
	glm::vec3 calculate_du(float x, float y, float z) const
	{
		return glm::normalize(glm::vec3(-max_phi*y, max_phi*x, 0));
	}
	glm::vec3 calculate_dv(float u, float v) const
	{
		return glm::normalize(glm::vec3(0, 0, max_z - min_z));
	}
	glm::vec3 calculate_dv(float x, float y, float z) const
	{
		return glm::normalize(glm::vec3(0, 0, max_z - min_z));
	}
	
	glm::vec3 calculate_normal(float u, float v) const
	{
		return glm::normalize(glm::cross(calculate_du(u, v), calculate_dv(u, v)));
	}
	
	glm::vec3 calculate_normal(float x, float y, float z) const
	{
		//gradient f(x,y,z) = 0, (df/dx,df/dy,dy/dz)
		glm::vec3 N = glm::vec3(2 * x, 2 * y, 0);
		return glm::normalize(N);
	}
private:
	//implicit equation x^2 + y^2 - r^2 = 0
	//uv is  u is angle, v is height. phi = u*maxphi, x=r*cos(phi), y=r*sin(phi), z=min_z + v(max_z-min_z)

	float r;
	float min_z, max_z;
	float max_phi; //in radians
};

class Disk : public Shape
{
	struct DiskIntersect
	{
		float t;
		glm::vec3 hitp;
		glm::vec3 ray_d;
		float phi;
	};
public:
	Disk(const std::string& _name, glm::mat4 rigidtransform, float height, float in_radius, float out_radius, float _phimax)
		: Shape(_name, rigidtransform)
	{ 
		inner_r = in_radius;
		outer_r = out_radius;
		h = height;
		phimax = glm::radians(_phimax);
	}
	~Disk() = default;

	float Area() const 
	{ 
		return phimax*.5f*(outer_r*outer_r - inner_r*inner_r); 
	}
	
	Bounds3 Bounds() const 
	{ 
		return TransformBounds(Bounds3(glm::vec3(-outer_r, -outer_r, h), glm::vec3(outer_r, outer_r, h)),ObjectToRender);
	}
	
	std::optional<LocalSurfaceInfo> Intersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		std::optional<DiskIntersect> isect = BasicIntersect(ray, tMax);
		if (!isect)
			return {};

		glm::vec3 p = isect.value().hitp;

		LocalSurfaceInfo info;
		info.tHit = isect.value().t;
		info.hitp = p;
		glm::vec2 uv = PostoUV(p.x, p.y, p.z);
		info.u = glm::clamp(uv.x, 0.0f, 1.0f);;
		info.v = glm::clamp(uv.y, 0.0f, 1.0f);
		info.du = calculate_du(p.x, p.y, p.z);
		info.dv = calculate_dv(p.x, p.y, p.z);
		info.n = calculate_normal(p.x, p.y, p.z);
		if (glm::dot(info.n, isect.value().ray_d) > 0) //the normal and the ray direction should be on opposite sides
			info.n = -info.n;
		info.wo = glm::vec3(0, 0, 1);

		//transform to world space
		info.Transform(ObjectToRender);

		return LocalSurfaceInfo{ info };
	}

	bool IntersectP(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		return BasicIntersect(ray, tMax).has_value();
	}

	std::optional<DiskIntersect> BasicIntersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		//transform ray from rendertoobject space
		glm::vec3 o = RenderToObject * glm::vec4(ray.o.x, ray.o.y, ray.o.z, 1);
		glm::vec3 d = RenderToObject * glm::vec4(ray.d.x, ray.d.y, ray.d.z, 0);

		//find intersecton with plane z=height
		float t0 = (h - o.z) / d.z;
		if (t0 <= 0 || t0 >= tMax)
			return {};

		if (d.z == 0)
			return {};

		//then see if distance from axis is within inner and outer radius and within the phimax
		glm::vec3 phit = o + t0 * d;
		float dist2 = phit.x * phit.x + phit.y * phit.y;
		if (dist2 > outer_r * outer_r || dist2 < inner_r * inner_r)
			return {};

		float phi = std::atan2(phit.y, phit.x);
		if (phi < 0) phi += 2 * std::numbers::pi;
		if (phi > phimax)
			return {};

		return DiskIntersect{t0, phit, glm::normalize(d), phi};
	}

	glm::vec3 UVtoPos(float u, float v) const
	{
		float phi = u * phimax;
		glm::vec3 pos(((1 - v) * outer_r + v * inner_r) * cos(phi), ((1 - v) * outer_r + v * inner_r) * sin(phi), h);

		return pos;
	}

	//assuming x,y,z is on the shape
	glm::vec2 PostoUV(float x, float y, float z) const
	{
		float phi = std::atan2(y, x);
		if (phi < 0) phi += 2 * std::numbers::pi;
		float u = phi / phimax;
		float v = (outer_r - std::sqrt(x * x + y * y)) / (outer_r - inner_r);

		return glm::vec2(u, v);
	}

	glm::vec3 calculate_du(float x, float y, float z) const
	{
		return glm::normalize(glm::vec3(-phimax * y, phimax * x, 0));
	}

	glm::vec3 calculate_dv(float x, float y, float z) const
	{
		return glm::normalize(glm::vec3(x, y, 0) * (inner_r - outer_r) / std::sqrt(x * x + y * y));
	}

	glm::vec3 calculate_normal(float u, float v) const
	{
		return glm::vec3(0, 0, 1);
	}

	glm::vec3 calculate_normal(float x, float y, float z) const
	{
		return glm::vec3(0, 0, 1);
	}
private:
	//no implicit form right now
	//disk parameterization, x= ((1-v)r + vri)cos(phi), y=((1-v)r+vri)sin(phi), z=h, phi=u*phimax

	float inner_r;
	float outer_r;
	float h;
	float phimax;
};

class TriangleSimple : public Shape
{
	struct TriangleIntersect
	{
		float t;
		glm::vec3 hitp;
		glm::vec3 ray_d;
		float B, Y;
	};

public:
	TriangleSimple(const std::string& _name, glm::mat4 rigidtransform, glm::vec3 _p1, glm::vec3 _p2, glm::vec3 _p3) : Shape(_name, rigidtransform)
	{ 
		p1 = _p1;
		p2 = _p2;
		p3 = _p3;
	}
	~TriangleSimple() = default;

	float Area() const 
	{ 
		return	0.5f * glm::length(glm::cross(p2 - p1, p3 - p1));
	}

	Bounds3 Bounds() const 
	{
		return TransformBounds(Bounds3(glm::vec3(glm::min(glm::min(p1.x, p2.x), p3.x),
					             glm::min(glm::min(p1.y, p2.y), p3.y), 
			                     glm::min(glm::min(p1.z, p2.z), p3.z)), 
						glm::vec3(glm::max(glm::max(p1.x, p2.x), p3.x), 
								  glm::max(glm::max(p1.y, p2.y), p3.y), 
								  glm::max(glm::max(p1.z, p2.z), p3.z))), ObjectToRender);
	}

	std::optional<LocalSurfaceInfo> Intersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		std::optional<TriangleIntersect> isect = BasicIntersect(ray, tMax);
		if (!isect)
			return {};

		glm::vec3 p = isect.value().hitp;

		LocalSurfaceInfo info;
		info.tHit = isect.value().t;
		info.hitp = p;
		//I can do barycentric coordinates
		//glm::vec2 uv = PostoUV(p.x, p.y, p.z);
		info.u = glm::clamp(isect.value().B,0.0f,1.0f);
		info.v = glm::clamp(isect.value().Y, 0.0f, 1.0f);
		info.du = calculate_du();
		info.dv = calculate_dv();
		info.n = calculate_normal();
		if (glm::dot(info.n, isect.value().ray_d) > 0) //the normal and the ray direction should be on opposite sides
			info.n = -info.n;
		info.wo = glm::vec3(0, 0, 1);

		//transform to world space
		info.Transform(ObjectToRender);

		return LocalSurfaceInfo{ info };
	}

	bool IntersectP(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		return BasicIntersect(ray, tMax).has_value();
	}

	//intersection 2 converts it to ray space, then you have a distance-type function and you see if the point lies on the same side of all 3 sides
	//of the triangle

	std::optional<TriangleIntersect> BasicIntersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		Timer timer1;
		timer1.Begin();
		//transform ray from rendertoobject space
		glm::vec3 orig = RenderToObject * glm::vec4(ray.o.x, ray.o.y, ray.o.z, 1);
		glm::vec3 dir = RenderToObject * glm::vec4(ray.d.x, ray.d.y, ray.d.z, 0);

		float a = p1.x - p2.x;
		float b = p1.y - p2.y;
		float c = p1.z - p2.z;
		float d = p1.x - p3.x;
		float e = p1.y - p3.y;
		float f = p1.z - p3.z;
		float g = dir.x;
		float h = dir.y;
		float i = dir.z;
		float j = p1.x - orig.x;
		float k = p1.y - orig.y;;
		float l = p1.z - orig.z;

		float M = a * (e * i - h * f) + b * (g * f - d * i) + c * (d * h - e * g);
		
		float t = -(f * (a * k - j * b) + e * (j * c - a * l) + d * (b * l - k * c)) / M;
		if (t < 0 || t >= tMax)
			return {};
		
		float Y = (i * (a * k - j * b) + h * (j * c - a * l) + g * (b * l - k * c)) / M;
		if (Y < 0 || Y > 1)
			return {};

		float B = (j * (e * i - h * f) + k * (g * f - d * i) + l * (d * h - e * g)) / M;
		if (B < 0 || B > 1- Y)
			return {};
		

		glm::vec3 hitp = orig + t * dir;
		//std::cout << "simple time: " << timer1.getTimeMicro() << std::endl;
		return TriangleIntersect{t, hitp, glm::normalize(dir), B, Y};
	}

	glm::vec3 UVtoPos(float u, float v) const
	{
		glm::vec3 pos = p1 + (p2 - p1) * u + (p3 - p1) * v;

		return pos;
	}

	//assuming x,y,z is on the shape
	//glm::vec2 PostoUV(float x, float y, float z) const
	//{
		//return glm::vec2(0, 0);
	//}

	glm::vec3 calculate_du() const
	{
		return glm::normalize(p2 - p1);
	}

	glm::vec3 calculate_dv() const
	{
		return glm::normalize(p3-p1);
	}

	glm::vec3 calculate_normal() const
	{
		//orientation??
		return glm::normalize(glm::cross(p3 - p1, p2 - p1));
	}

private:
	//parameterization of triangle we have B and Y, position = p1 + (p2-p1)*B + (p3-p1)*Y, B>0 and Y>0 and B + Y < 1
	//x = p1.x + (p2.x-p1.x)*B + (p3.x-p1.x)*Y, y and z are same

	glm::vec3 p1;
	glm::vec3 p2;
	glm::vec3 p3;
};

namespace Hitdata {
	static int triangle_intersect_count = 0;
}

class Triangle : public Shape
{
public:
	//do we grab this info from the data, or just calculate it outselves?
	struct vertex_available
	{
		bool texcoords = true;
		bool normals = true;
		bool tangents = true;
		bool bitangents = true;
		bool precomputed_worldtransform = false;
	};

	struct TriangleIntersect
	{
		//bary coordinates
		float b0;
		float b1;
		float b2;
		float t;
		glm::vec3 rayd;
	};

public:
	//**these are given in object space, some of the algorithms I think assume this is in world space
	//calculations are done in world space, ray and points expected to be in world space
	Triangle(const std::string& _name, glm::mat4 rigidtransform, std::string _model_name, int _mesh_id, int _tri_id, 
		vertex_available avail_info = vertex_available()) : Shape(_name, rigidtransform)
	{
		model_name = _model_name;
		mesh_id = _mesh_id;
		tri_id = _tri_id;
		available_info = avail_info;
	}
	~Triangle() = default;

	float Area() const
	{
		if (modelcache.contains(model_name) == true)
		{
			MeshCache::Mesh& mesh = modelcache[model_name].meshes[mesh_id];
			glm::vec3 p0 = mesh.positions[mesh.indices[3 * tri_id]];
			glm::vec3 p1 = mesh.positions[mesh.indices[3 * tri_id + 1]];
			glm::vec3 p2 = mesh.positions[mesh.indices[3 * tri_id + 2]];
			
			return	0.5f * glm::length(glm::cross(p1 - p0, p2 - p0));
		}
		return 0.0f;
	}

	Bounds3 Bounds() const
	{
		if (modelcache.contains(model_name) == true)
		{
			MeshCache::Mesh& mesh = modelcache[model_name].meshes[mesh_id];
			glm::vec3 p0 = mesh.positions[mesh.indices[3 * tri_id]];
			glm::vec3 p1 = mesh.positions[mesh.indices[3 * tri_id + 1]];
			glm::vec3 p2 = mesh.positions[mesh.indices[3 * tri_id + 2]];

			return TransformBounds(Bounds3(glm::vec3(glm::min(glm::min(p0.x, p1.x), p2.x),
				glm::min(glm::min(p0.y, p1.y), p2.y),
				glm::min(glm::min(p0.z, p1.z), p2.z)),
				glm::vec3(glm::max(glm::max(p0.x, p1.x), p2.x),
					glm::max(glm::max(p0.y, p1.y), p2.y),
					glm::max(glm::max(p0.z, p1.z), p2.z))), ObjectToRender);
		}
		return Bounds3(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0));
	}

	std::optional<LocalSurfaceInfo> CalculateLocalSurface(TriangleIntersect isect) const
	{
		//everythings basically interpolating barycentric coordinates based on values in the 3 vertices
		if (modelcache.contains(model_name) == false)
		{
			std::cout << "not found in model cache triintersect\n";
			return {};
		}
		MeshCache::Mesh& mesh = modelcache[model_name].meshes[mesh_id];
		glm::vec3 p0 = mesh.positions[mesh.indices[3 * tri_id]];
		glm::vec3 p1 = mesh.positions[mesh.indices[3 * tri_id + 1]];
		glm::vec3 p2 = mesh.positions[mesh.indices[3 * tri_id + 2]];

		glm::vec3 p0_w = ObjectToRender * glm::vec4(p0.x, p0.y, p0.z, 1);
		glm::vec3 p1_w = ObjectToRender * glm::vec4(p1.x, p1.y, p1.z, 1);
		glm::vec3 p2_w = ObjectToRender * glm::vec4(p2.x, p2.y, p2.z, 1);

		//for now (u,v) p0(0,0), p1(1,0), p2(0,1)
		std::array<glm::vec2, 3> uv_values = { glm::vec2(0, 0), glm::vec2(1, 0), glm::vec2(0, 1) };

		glm::vec2 duv02 = uv_values[0] - uv_values[2], duv12 = uv_values[1] - uv_values[2];
		glm::vec3 dp02 = p0_w - p2_w, dp12 = p1_w - p2_w;
		float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];//DifferenceOfProducts(duv02[0], duv12[1], duv02[1], duv12[0]);

		glm::vec3 dpdu, dpdv;
		bool degenerateUV = std::abs(determinant) < 1e-9f;
		if (!degenerateUV) {
			// Compute triangle $\dpdu$ and $\dpdv$ via matrix inversion
			float invdet = 1 / determinant;
			dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;// DifferenceOfProducts(duv12[1], dp02, duv02[1], dp12)* invdet;
			dpdv = (duv02[0] * dp12 - duv12[0] * dp02) * invdet;// DifferenceOfProducts(duv02[0], dp12, duv12[0], dp02)* invdet;
		}
		// Handle degenerate triangle $(u,v)$ parameterization or partial derivatives
		if (degenerateUV || std::pow(glm::length(glm::cross(dpdu, dpdv)), 2) == 0) {
			glm::vec3 ng = glm::cross(p2_w - p0_w, p1_w - p0_w);
			if (std::pow(glm::length(ng), 2) == 0) {
				ng = glm::vec3(glm::cross(glm::dvec3(p2_w - p0_w), glm::dvec3(p1_w - p0_w)));
				if (glm::length(ng) == 0)
					std::cout << "triangle interaction ng == 0";
				//CHECK_NE(LengthSquared(ng), 0);
			}
			ng = glm::normalize(ng);
			//CoordinateSystem(Normalize(ng), &dpdu, &dpdv);
			float sign = std::copysign(float(1), ng.z);
			float a = -1 / (sign + ng.z);
			float b = ng.x * ng.y * a;
			dpdu = glm::vec3(1 + sign * std::pow(ng.x, 2) * a, sign * b, -sign * ng.x);
			dpdv = glm::vec3(b, sign + std::pow(ng.y, 2) * a, -ng.y);
		}

		LocalSurfaceInfo info;
		//pos
		info.hitp = p0_w * isect.b0 + p1_w * isect.b1 + p2_w * isect.b2;

		//uv
		glm::vec2 tc1 = mesh.texcoords[mesh.indices[3 * tri_id]];
		glm::vec2 tc2 = mesh.texcoords[mesh.indices[3 * tri_id + 1]];
		glm::vec2 tc3 = mesh.texcoords[mesh.indices[3 * tri_id + 2]];
		glm::vec2 hit_uv;
		if (available_info.texcoords)
			hit_uv = tc1 * isect.b0 + tc2 * isect.b1 + tc3 * isect.b2;
		else
			hit_uv = uv_values[0] * isect.b0 + uv_values[1] * isect.b1 + uv_values[2] * isect.b2;
		info.u = glm::clamp(hit_uv.x, 0.0f, 1.0f);
		info.v = glm::clamp(hit_uv.y, 0.0f, 1.0f);

		//tangents
		glm::vec3 tan1 = mesh.tangents[mesh.indices[3 * tri_id]];
		glm::vec3 tan2 = mesh.tangents[mesh.indices[3 * tri_id + 1]];
		glm::vec3 tan3 = mesh.tangents[mesh.indices[3 * tri_id + 2]];
		if (available_info.tangents)
			info.du = glm::normalize(tan1 * isect.b0 + tan2 * isect.b1 + tan3 * isect.b2);
		else
			info.du = dpdu;

		//bitangents
		glm::vec3 bitan1 = mesh.bitangents[mesh.indices[3 * tri_id]];
		glm::vec3 bitan2 = mesh.bitangents[mesh.indices[3 * tri_id + 1]];
		glm::vec3 bitan3 = mesh.bitangents[mesh.indices[3 * tri_id + 2]];
		if (available_info.bitangents)
			info.dv = glm::normalize(bitan1 * isect.b0 + bitan2 * isect.b1 + bitan3 * isect.b2);
		else
			info.dv = dpdv;

		//normal
		glm::vec3 n1 = mesh.normals[mesh.indices[3 * tri_id]];
		glm::vec3 n2 = mesh.normals[mesh.indices[3 * tri_id + 1]];
		glm::vec3 n3 = mesh.normals[mesh.indices[3 * tri_id + 2]];
		if (available_info.normals)
			info.n = glm::normalize(n1 * isect.b0 + n2 * isect.b1 + n3 * isect.b2);
		else
			info.n = glm::normalize(glm::cross(dp02, dp12));
		if (glm::dot(info.n, isect.rayd) > 0) //the normal and the ray direction should be on opposite sides
			info.n = -info.n;

		//wo
		info.wo = glm::vec3(0, 0, 1);

		//already in world space

		return LocalSurfaceInfo{ info };
	}

	std::optional<LocalSurfaceInfo> Intersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		std::optional<TriangleIntersect> opt_isect = BasicIntersect(ray, tMax);
		if (!opt_isect)
			return {};
		
		return CalculateLocalSurface(opt_isect.value());
	}

	bool IntersectP(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		return BasicIntersect(ray, tMax).has_value();
	}

	//static int triangle_intersect_count;
	//pbr implementation
	std::optional<TriangleIntersect> BasicIntersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		if (modelcache.contains(model_name) == false)
		{
			std::cout << "tri basic intersect cant find model\n";
			return {};
		}
		Hitdata::triangle_intersect_count++;
		Timer timer1;
		timer1.Begin();
		MeshCache::Mesh& mesh = modelcache[model_name].meshes[mesh_id];
		glm::vec3 p0 = mesh.positions[mesh.indices[3 * tri_id]];
		glm::vec3 p1 = mesh.positions[mesh.indices[3 * tri_id + 1]];
		glm::vec3 p2 = mesh.positions[mesh.indices[3 * tri_id + 2]];

		glm::vec3 p0_w, p1_w, p2_w;
		if (!available_info.precomputed_worldtransform)
		{
			p0_w = ObjectToRender * glm::vec4(p0.x, p0.y, p0.z, 1);
			p1_w = ObjectToRender * glm::vec4(p1.x, p1.y, p1.z, 1);
			p2_w = ObjectToRender * glm::vec4(p2.x, p2.y, p2.z, 1);
		}
		else
		{
			p0_w = p0;
			p1_w = p1;
			p2_w = p2;
		}

		 // Return no intersection if triangle is degenerate
		if (std::pow(glm::length(glm::cross(p2_w - p0_w, p1_w - p0_w)), 2) == 0) {
			//std::cout << "degen\n";
			return {};
		}
		
		// Transform triangle vertices to ray coordinate space
		// Translate vertices based on ray origin
		glm::vec3 p0t = p0_w - glm::vec3(ray.o);
		glm::vec3 p1t = p1_w - glm::vec3(ray.o);
		glm::vec3 p2t = p2_w - glm::vec3(ray.o);

		// Permute components of triangle vertices and ray direction
		glm::vec3 positive_d = glm::vec3(std::abs(ray.d.x), std::abs(ray.d.y), std::abs(ray.d.z));

		int kz = pbrt::MaxComponentIndex(positive_d); //max_index; //MaxComponentIndex(Abs(ray.d));
		int kx = kz + 1;
		if (kx == 3)
			kx = 0;
		int ky = kx + 1;
		if (ky == 3)
			ky = 0;
		glm::vec3 d = glm::vec3(ray.d[kx], ray.d[ky], ray.d[kz]); //Permute(ray.d, { kx, ky, kz });
		p0t = glm::vec3(p0t[kx], p0t[ky], p0t[kz]); //Permute(p0t, { kx, ky, kz });
		p1t = glm::vec3(p1t[kx], p1t[ky], p1t[kz]); // Permute(p1t, { kx, ky, kz });
		p2t = glm::vec3(p2t[kx], p2t[ky], p2t[kz]); //Permute(p2t, { kx, ky, kz });

		// Apply shear transformation to translated vertex positions
		float Sx = -d.x / d.z;
		float Sy = -d.y / d.z;
		float Sz = 1 / d.z;
		p0t.x += Sx * p0t.z;
		p0t.y += Sy * p0t.z;
		p1t.x += Sx * p1t.z;
		p1t.y += Sy * p1t.z;
		p2t.x += Sx * p2t.z;
		p2t.y += Sy * p2t.z;

		// Compute edge function coefficients _e0_, _e1_, and _e2_
		float e0 = pbrt::DifferenceOfProducts(p1t.x, p2t.y, p1t.y, p2t.x); //(p1t.x*p2t.y - p1t.y*p2t.x) + (-p1t.y*p2t.x + p1t.y*p2t.x); // DifferenceOfProducts(p1t.x, p2t.y, p1t.y, p2t.x);
		float e1 = pbrt::DifferenceOfProducts(p2t.x, p0t.y, p2t.y, p0t.x); //(p2t.x * p0t.y - p2t.y * p0t.x) + (-p2t.y * p0t.x + p2t.y * p0t.x); //DifferenceOfProducts(p2t.x, p0t.y, p2t.y, p0t.x);
		float e2 = pbrt::DifferenceOfProducts(p0t.x, p1t.y, p0t.y, p1t.x); //(p0t.x * p1t.y - p0t.y * p1t.x) + (-p0t.y * p1t.x + p0t.y * p1t.x); //DifferenceOfProducts(p0t.x, p1t.y, p0t.y, p1t.x);

		// Fall back to double-precision test at triangle edges
		if (sizeof(float) == sizeof(float) && (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
			double p2txp1ty = (double)p2t.x * (double)p1t.y;
			double p2typ1tx = (double)p2t.y * (double)p1t.x;
			e0 = (float)(p2typ1tx - p2txp1ty);
			double p0txp2ty = (double)p0t.x * (double)p2t.y;
			double p0typ2tx = (double)p0t.y * (double)p2t.x;
			e1 = (float)(p0typ2tx - p0txp2ty);
			double p1txp0ty = (double)p1t.x * (double)p0t.y;
			double p1typ0tx = (double)p1t.y * (double)p0t.x;
			e2 = (float)(p1typ0tx - p1txp0ty);
		}

		// Perform triangle edge and determinant tests
		if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0)) {
			//std::cout << "e's\n";
			return {};
		}
		float det = e0 + e1 + e2;
		if (det == 0) {
			//std::cout << "det 0\n";
			return {};
		}

		// Compute scaled hit distance to triangle and test against ray $t$ range
		p0t.z *= Sz;
		p1t.z *= Sz;
		p2t.z *= Sz;
		float tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
		if (det < 0 && (tScaled >= 0 || tScaled < tMax * det)) {
			//std::cout << "3\n";
			return {};
		}

		else if (det > 0 && (tScaled <= 0 || tScaled > tMax * det)) {
			//std::cout << "4\n";
			return {};
		}

		// Compute barycentric coordinates and $t$ value for triangle intersection
		float invDet = 1 / det;
		float b0 = e0 * invDet, b1 = e1 * invDet, b2 = e2 * invDet;
		float t = tScaled * invDet;
		if (std::isnan(t))
		{
			//std::cout << tScaled << " " << invDet << std::endl;
			return {};
			//std::cout << "triangle intersection isnan(t)";
		}
		//DCHECK(!IsNaN(t));

		// Ensure that computed triangle $t$ is conservatively greater than zero
		// Compute $\delta_z$ term for triangle $t$ error bounds
		glm::vec3 positive_p = glm::vec3(std::abs(p0t.z), std::abs(p1t.z), std::abs(p2t.z));
		//max_index = 0; (positive_p.y > positive_p.x) ? max_index = 1 : max_index = 0; (positive_p.z > positive_p[max_index]) ? max_index = 2 : max_index = 0;
		float maxZt = pbrt::MaxComponentValue(positive_p);// positive_p[max_index];// MaxComponentValue(Abs(glm::vec3(p0t.z, p1t.z, p2t.z)));
		float deltaZ = pbrt::gamma(3) * maxZt;

		// Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
		positive_p = glm::vec3(std::abs(p0t.x), std::abs(p1t.x), std::abs(p2t.x));
		//max_index = 0; (positive_p.y > positive_p.x) ? max_index = 1 : max_index = 0; (positive_p.z > positive_p[max_index]) ? max_index = 2 : max_index = 0;
		float maxXt = pbrt::MaxComponentValue(positive_p);// MaxComponentValue(Abs(Vector3f(p0t.x, p1t.x, p2t.x)));
		
		positive_p = glm::vec3(std::abs(p0t.y), std::abs(p1t.y), std::abs(p2t.y));
		//max_index = 0; (positive_p.y > positive_p.x) ? max_index = 1 : max_index = 0; (positive_p.z > positive_p[max_index]) ? max_index = 2 : max_index = 0;
		float maxYt = pbrt::MaxComponentValue(positive_p); // MaxComponentValue(Abs(Vector3f(p0t.y, p1t.y, p2t.y)));
		float deltaX = pbrt::gamma(5) * (maxXt + maxZt);
		float deltaY = pbrt::gamma(5) * (maxYt + maxZt);

		// Compute $\delta_e$ term for triangle $t$ error bounds
		float deltaE = 2 * (pbrt::gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

		// Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
		positive_p = glm::vec3(std::abs(e0), std::abs(e1), std::abs(e2));
		//max_index = 0; (positive_p.y > positive_p.x) ? max_index = 1 : max_index = 0; (positive_p.z > positive_p[max_index]) ? max_index = 2 : max_index = 0;
		float maxE = pbrt::MaxComponentValue(positive_p); ;// positive_p[max_index];// MaxComponentValue(Abs(Vector3f(e0, e1, e2)));
		float deltaT =
			3 * (pbrt::gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) * std::abs(invDet);
		if (t <= deltaT) {
			// std::cout << "5 " << t << " " << deltaT << " " << invDet << " z: " << p0t.z << " " << ray.o.z << std::endl;
			//std::cout << ray.o.x << " " << ray.o.y << " " << ray.o.z << " dir: " << ray.d.x << " " << ray.d.y << " " << ray.d.z 
				//<<" p0t: "<<p0t.x << " " << p0t.y << " " << p0t.z << std::endl;
			return {};
		}
		//std::cout << "time: " << timer1.getTimeMicro() << std::endl;
		// Return _TriangleIntersection_ for intersection
		return TriangleIntersect{ b0, b1, b2, t, glm::normalize(ray.d) };
	}

private:
	std::string model_name;
	int mesh_id;
	int tri_id;
	std::unordered_map<std::string, MeshCache::Model>& modelcache = MeshCache::modelCache;
	
	vertex_available available_info;
};


class TriModel : public Shape
{
	struct TriModelIntersect
	{
		Triangle::TriangleIntersect tri_isect;
		int mesh_id;
		int tri_id;
	};

public:
	TriModel(const std::string& _name, glm::mat4 rigidtransform, std::string _model_name, bool back_facing_cull, bool precomp_worldspace_model, Triangle::vertex_available avail)
		: precomputed_bounds(glm::vec3(0,0,0),glm::vec3(0,0,0)), mesh_has_precomputed_worldposition(precomp_worldspace_model),
		  model_name(_model_name), enable_cull_back_face(back_facing_cull), avail_info(avail), Shape(_name, rigidtransform)
	{
		if (!modelcache.contains(model_name))
			std::cout << "trimodel model doesnt exist\n";	
		MeshCache::Model& model = modelcache[model_name];

		//generate precomputed bounds,  find min/max (x,y,z) out of every vertex
		glm::vec3 min_p = glm::vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
		glm::vec3 max_p = glm::vec3(std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
		for (int mesh_id = 0; mesh_id < model.meshes.size(); mesh_id++)
		{
			MeshCache::Mesh& mesh = model.meshes[mesh_id];
			for (int i = 0; i < mesh.positions.size(); i++)
			{
				glm::vec3 pos = mesh.positions[i];
				min_p.x = std::min(min_p.x, pos.x); min_p.y = std::min(min_p.y, pos.y); min_p.z = std::min(min_p.z, pos.z);
				max_p.x = std::max(max_p.x, pos.x); max_p.y = std::max(max_p.y, pos.y); max_p.z = std::max(max_p.z, pos.z);
			}
		}
		precomputed_bounds.pmin = min_p;
		precomputed_bounds.pmax = max_p;
		//transform to world?

		//ComputeBackFace(glm::vec3(0, 0, 1));
		
		//generate triangles
		triangles.resize(model.meshes.size());
		for (int mesh_id = 0; mesh_id < model.meshes.size(); mesh_id++)
		{
			MeshCache::Mesh& mesh = model.meshes[mesh_id];
			triangles[mesh_id].reserve(mesh.indices.size() / 3);
			for (int tri_id = 0; tri_id < mesh.indices.size() / 3; tri_id++)
			{
				triangles[mesh_id].push_back(Triangle("tri", rigidtransform, model_name, mesh_id, tri_id, avail_info));
			}
		}
	}
	~TriModel() = default;

	void SetRigidTransform(glm::mat4 rigidtransform) override
	{
		//go through every triangle and call
		for (int i = 0; i < triangles.size(); i++)
		{
			for (int j = 0; j < triangles[i].size(); j++)
			{
				triangles[i][j].SetRigidTransform(rigidtransform);
			}
		}
		//call for oursevles (bounding box)
		Shape::SetRigidTransform(rigidtransform);
	}

	//if you move look position or move the model well you have to recalculate this
	void EnableBackface(bool val) { enable_cull_back_face = val; }
	void ComputeBackFace(glm::vec3 world_pos_look, bool enable_backface)
	{
		enable_cull_back_face = enable_backface;
		if (enable_cull_back_face)
		{
			MeshCache::Model& model = modelcache[model_name];

			glm::vec3 look_dir = glm::normalize(world_pos_look);

			back_facing.clear();
			back_facing.resize(model.meshes.size());
			int counter = 0;
			for (int mesh_id = 0; mesh_id < model.meshes.size(); mesh_id++)
			{
				MeshCache::Mesh& mesh = model.meshes[mesh_id];
				back_facing[mesh_id].reserve(mesh.indices.size() / 3);
				for (int tri_id = 0; tri_id < mesh.indices.size() / 3; tri_id++)
				{
					//take the average of the face (could also just compute it with pos)
					glm::vec3 n1 = mesh.normals[mesh.indices[3 * tri_id]];
					glm::vec3 n2 = mesh.normals[mesh.indices[3 * tri_id + 1]];
					glm::vec3 n3 = mesh.normals[mesh.indices[3 * tri_id + 2]];
					glm::vec3 N = glm::normalize((n1 + n2 + n3) / 3.0f);
					if (!mesh_has_precomputed_worldposition)
					{
						glm::mat3 normal_transform = glm::transpose(glm::inverse(ObjectToRender));
						N = glm::normalize(normal_transform * N);
					}

					//world pos
					if (glm::dot(look_dir, N) > 0) {
						back_facing[mesh_id].push_back(true);
						counter++;
					}
					else {
						back_facing[mesh_id].push_back(false);
					}
				}
			}
			//std::cout << counter << std::endl;
		}	
	}

	float Area() const
	{
		//estimate with bounds
		return  (precomputed_bounds.pmax.x - precomputed_bounds.pmin.x) 
			  * (precomputed_bounds.pmax.y - precomputed_bounds.pmin.y)
			  * (precomputed_bounds.pmax.z - precomputed_bounds.pmin.z);
	}

	Bounds3 Bounds() const
	{
		//precompute it in model space and then transform to world whenever, unless points are precomputed in world
		if(!mesh_has_precomputed_worldposition)
			return TransformBounds(precomputed_bounds, ObjectToRender);
		else
			return precomputed_bounds;
	}

	std::optional<LocalSurfaceInfo> Intersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		std::optional<TriModelIntersect> opt_isect = BasicIntersect(ray, tMax);
		if (!opt_isect)
			return {};
		TriModelIntersect isect = opt_isect.value();

		return triangles[isect.mesh_id][isect.tri_id].CalculateLocalSurface(isect.tri_isect);
	}

	bool IntersectP(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		return BasicIntersect(ray, tMax).has_value();
	}

	std::optional<TriModelIntersect> BasicIntersect(const Ray& ray, float tMax = std::numeric_limits<float>::max()) const
	{
		//go through every triangle keep the closest t value, and thats the one we give back
		//first check the bounding intersect, then check cullface if we have that optimization
		
		//check intersection with precomputed bounding box
		Bounds3 bounding_box = precomputed_bounds;
		if (!mesh_has_precomputed_worldposition)
			bounding_box = TransformBounds(precomputed_bounds, ObjectToRender);
		
		if (!bounding_box.IntersectP(ray, tMax))
		{
			return {};
		}
		
		//go through every triangle and get closest hit_t
		float t_min = std::numeric_limits<float>::max();
		int id_mesh = 0;
		int id_tri = 0;
		std::optional<Triangle::TriangleIntersect> current_isect;
		MeshCache::Model& model = modelcache[model_name];
		int counter = 0;
		for (int mesh_id = 0; mesh_id < triangles.size(); mesh_id++)
		{
			const std::vector<Triangle>& tris = triangles[mesh_id];
			for (int tri_id = 0; tri_id < tris.size(); tri_id++)
			{
				//early leave if its backfacing
				if (enable_cull_back_face && !back_facing.empty() && back_facing[mesh_id][tri_id]) {
					//std::cout << "culled out\n";
					continue;
				}
				if (triangles[mesh_id][tri_id].BasicIntersect(ray, tMax).has_value())
				{
					counter++;
					std::optional<Triangle::TriangleIntersect> isect = triangles[mesh_id][tri_id].BasicIntersect(ray, tMax);
					if (isect.has_value() && (isect.value().t < t_min))
					{
						id_mesh = mesh_id;
						id_tri = tri_id;
						t_min = isect.value().t;
						current_isect = isect;
					}
				}
				//optimization if I want, will break tho at cracks
				//if (counter >= 1)
					//break;
			}
		}
		//std::cout << "hit: " << counter << std::endl;
		if (t_min == std::numeric_limits<float>::max())
			return {};

		//everything should be world space, triangles work in world space

		return TriModelIntersect{ current_isect.value(), id_mesh, id_tri};

	}

	std::string getModelName() { return model_name; }
private:
	friend class Octtree_Model;

	//just groups triangles together as one model loaded, has some higher level functionality
	//bounding box, backwards facing removal, precomputed world_position
	
	std::string model_name;
	Triangle::vertex_available avail_info;
	std::vector<std::vector<Triangle>> triangles; //for every mesh, it holds every triangle
	Bounds3 precomputed_bounds; 

	bool enable_cull_back_face;
	std::vector<std::vector<bool>> back_facing; //true or false if its backfacing for quick intersection exit
	
	bool mesh_has_precomputed_worldposition; //is this model precomputed world position? (has to be done at the meshcache level)

	std::unordered_map<std::string, MeshCache::Model>& modelcache = MeshCache::modelCache;
};


//triangle mesh

//bilinear patch

//curve





