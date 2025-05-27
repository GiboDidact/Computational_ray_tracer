/********************************************************/

/* AABB-triangle overlap test code                      */

/* by Tomas Akenine-Möller                              */

/* Function: int triBoxOverlap(float boxcenter[3],      */

/*          float boxhalfsize[3],float triverts[3][3]); */

/* History:                                             */

/*   2001-03-05: released the code in its first version */

/*   2001-06-18: changed the order of the tests, faster */

/*                                                      */

/* Acknowledgement: Many thanks to Pierre Terdiman for  */

/* suggestions and discussions on how to optimize code. */

/* Thanks to David Hunt for finding a ">="-bug!         */

/********************************************************/

/*
Copyright 2020 Tomas Akenine-Möller

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and
to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "../pch.h"
#include <math.h>
#include <stdio.h>


/*
#define X 0

#define Y 1

#define Z 2



#define CROSS(dest,v1,v2) \

dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \

dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \

dest[2] = v1[0] * v2[1] - v1[1] * v2[0];



#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])



#define SUB(dest,v1,v2) \

dest[0] = v1[0] - v2[0]; \

dest[1] = v1[1] - v2[1]; \

dest[2] = v1[2] - v2[2];



#define FINDMINMAX(x0,x1,x2,min,max) \

min = max = x0;   \

if (x1 < min) min = x1; \

	if (x1 > max) max = x1; \

		if (x2 < min) min = x2; \

			if (x2 > max) max = x2;



#define AXISTEST_X01(a, b, fa, fb)			   \

p0 = a * v0[Y] - b * v0[Z];			       	   \

p2 = a * v2[Y] - b * v2[Z];			       	   \

if (p0 < p2) { min = p0; max = p2; }
else { min = p2; max = p0; } \

rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \

if (min > rad || max < -rad) return 0;



#define AXISTEST_X2(a, b, fa, fb)			   \

p0 = a * v0[Y] - b * v0[Z];			           \

p1 = a * v1[Y] - b * v1[Z];			       	   \

if (p0 < p1) { min = p0; max = p1; }
else { min = p1; max = p0; } \

rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \

if (min > rad || max < -rad) return 0;


#define AXISTEST_Y02(a, b, fa, fb)			   \

p0 = -a * v0[X] + b * v0[Z];		      	   \

p2 = -a * v2[X] + b * v2[Z];	       	       	   \

if (p0 < p2) { min = p0; max = p2; }
else { min = p2; max = p0; } \

rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \

if (min > rad || max < -rad) return 0;



#define AXISTEST_Y1(a, b, fa, fb)			   \

p0 = -a * v0[X] + b * v0[Z];		      	   \

p1 = -a * v1[X] + b * v1[Z];	     	       	   \

if (p0 < p1) { min = p0; max = p1; }
else { min = p1; max = p0; } \

rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \

if (min > rad || max < -rad) return 0;



#define AXISTEST_Z12(a, b, fa, fb)			   \

p1 = a * v1[X] - b * v1[Y];			           \

p2 = a * v2[X] - b * v2[Y];			       	   \

if (p2 < p1) { min = p2; max = p1; }
else { min = p1; max = p2; } \

rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \

if (min > rad || max < -rad) return 0;



#define AXISTEST_Z0(a, b, fa, fb)			   \

p0 = a * v0[X] - b * v0[Y];				   \

p1 = a * v1[X] - b * v1[Y];			           \

if (p0 < p1) { min = p0; max = p1; }
else { min = p1; max = p0; } \

rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \

if (min > rad || max < -rad) return 0;
*/

namespace Moller {

	void FindMinMax(float x0, float x1, float x2, float& min, float& max)
	{
		min = max = x0;
		if (x1 < min) min = x1;
		if (x1 > max) max = x1;
		if (x2 < min) min = x2;
		if (x2 > max) max = x2;
	}

	int planeBoxOverlap(glm::vec3 normal, glm::vec3 vert, glm::vec3 maxbox)	// -NJMP-

	{
		int q;

		glm::vec3 vmin, vmax;
		float v;

		for (q = 0; q <= 2; q++)
		{
			v = vert[q];					// -NJMP-
			if (normal[q] > 0.0f)
			{
				vmin[q] = -maxbox[q] - v;	// -NJMP-
				vmax[q] = maxbox[q] - v;	// -NJMP-
			}
			else
			{
				vmin[q] = maxbox[q] - v;	// -NJMP-
				vmax[q] = -maxbox[q] - v;	// -NJMP-
			}
		}

		if (glm::dot(normal, vmin) > 0.0f) return 0;	// -NJMP-

		if (glm::dot(normal, vmax) >= 0.0f) return 1;	// -NJMP-



		return 0;

	}

	int triBoxOverlap(glm::vec3 boxcenter, glm::vec3 boxhalfsize, std::array<glm::vec3, 3> triverts)
	{
		/*    use separating axis theorem to test overlap between triangle and box */

		/*    need to test for overlap in these directions: */

		/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */

		/*       we do not even need to test these) */

		/*    2) normal of the triangle */

		/*    3) crossproduct(edge from tri, {x,y,z}-directin) */

		/*       this gives 3x3=9 more tests */

		glm::vec3 v0, v1, v2;

		//   float axis[3];

		float min, max, p0, p1, p2, rad, fex, fey, fez;		// -NJMP- "d" local variable removed

		glm::vec3 normal, e0, e1, e2;



		/* This is the fastest branch on Sun */

		/* move everything so that the boxcenter is in (0,0,0) */
		v0 = triverts[0] - boxcenter;
		v1 = triverts[1] - boxcenter;
		v2 = triverts[2] - boxcenter;

		//SUB(v0, triverts[0], boxcenter);
		//SUB(v1, triverts[1], boxcenter);
		//SUB(v2, triverts[2], boxcenter);



		/* compute triangle edges */
		e0 = v1 - v0;
		e1 = v2 - v1;
		e2 = v0 - v2;

		//SUB(e0, v1, v0);      /* tri edge 0 */
		//SUB(e1, v2, v1);      /* tri edge 1 */
		//SUB(e2, v0, v2);      /* tri edge 2 */



		/* Bullet 3:  */

		auto AxisTest_X01 = [&](float a, float b, float fa, float fb) -> bool {
			p0 = a * v0.y - b * v0.z;
			p2 = a * v2.y - b * v2.z;

			if (p0 < p2) { min = p0; max = p2; }
			else { min = p2; max = p0; }

			rad = fa * boxhalfsize.y + fb * boxhalfsize.z;
			if (min > rad || max < -rad) return false;

			return true;
		};

		auto AxisTest_X2 = [&](float a, float b, float fa, float fb) -> bool {
			p0 = a * v0.y - b * v0.z;
			p1 = a * v1.y - b * v1.z;

			if (p0 < p1) { min = p0; max = p1; }
			else { min = p1; max = p0; }

			rad = fa * boxhalfsize.y + fb * boxhalfsize.z;
			if (min > rad || max < -rad) return false;

			return true;
		};

		auto AxisTest_Y02 = [&](float a, float b, float fa, float fb) -> bool {
			p0 = -a * v0.x + b * v0.z;

			p2 = -a * v2.x + b * v2.z;

			if (p0 < p2) { min = p0; max = p2; }
			else { min = p2; max = p0; }

			rad = fa * boxhalfsize.x + fb * boxhalfsize.z;
			if (min > rad || max < -rad) return false;

			return true;
		};

		auto AxisTest_Y1 = [&](float a, float b, float fa, float fb) -> bool {
			p0 = -a * v0.x + b * v0.z;
			p1 = -a * v1.x + b * v1.z;

			if (p0 < p1) { min = p0; max = p1; }
			else { min = p1; max = p0; }

			rad = fa * boxhalfsize.x + fb * boxhalfsize.z;
			if (min > rad || max < -rad) return false;

			return true;
		};

		auto AxisTest_Z0 = [&](float a, float b, float fa, float fb) -> bool {
			p0 = a * v0.x - b * v0.y;
			p1 = a * v1.x - b * v1.y;

			if (p0 < p1) { min = p0; max = p1; }
			else { min = p1; max = p0; }

			rad = fa * boxhalfsize.x + fb * boxhalfsize.y;
			if (min > rad || max < -rad) return true;

			return true;
		};

		auto AxisTest_Z12 = [&](float a, float b, float fa, float fb) -> bool {
			p1 = a * v1.x - b * v1.y;
			p2 = a * v2.x - b * v2.y;

			if (p2 < p1) { min = p2; max = p1; }
			else { min = p1; max = p2; }

			rad = fa * boxhalfsize.x + fb * boxhalfsize.y;
			if (min > rad || max < -rad) return false;

			return true;
		};


		/*  test the 9 tests first (this was faster) */

		fex = fabsf(e0.x);

		fey = fabsf(e0.y);

		fez = fabsf(e0.z);

		if (!AxisTest_X01(e0.z, e0.y, fez, fey))
			return false;

		if (!AxisTest_Y02(e0.z, e0.x, fez, fex))
			return false;

		if (!AxisTest_Z12(e0.y, e0.x, fey, fex))
			return false;

		//AXISTEST_X01(e0[Z], e0[Y], fez, fey);
		//AXISTEST_Y02(e0[Z], e0[X], fez, fex);
		//AXISTEST_Z12(e0[Y], e0[X], fey, fex);



		fex = fabsf(e1.x);

		fey = fabsf(e1.y);

		fez = fabsf(e1.z);

		if (!AxisTest_X01(e1.z, e1.y, fez, fey))
			return false;

		if (!AxisTest_Y02(e1.z, e1.x, fez, fex))
			return false;

		if (!AxisTest_Z0(e1.y, e1.x, fey, fex))
			return false;

		//AXISTEST_X01(e1[Z], e1[Y], fez, fey);
		//AXISTEST_Y02(e1[Z], e1[X], fez, fex);
		//AXISTEST_Z0(e1[Y], e1[X], fey, fex);



		fex = fabsf(e2.x);

		fey = fabsf(e2.y);

		fez = fabsf(e2.z);

		if (!AxisTest_X2(e2.z, e2.y, fez, fey))
			return false;

		if (!AxisTest_Y1(e2.z, e2.x, fez, fex))
			return false;

		if (!AxisTest_Z12(e2.y, e2.x, fey, fex))
			return false;

		//AXISTEST_X2(e2[Z], e2[Y], fez, fey);
		//AXISTEST_Y1(e2[Z], e2[X], fez, fex);
		//AXISTEST_Z12(e2[Y], e2[X], fey, fex);

		/* Bullet 1: */

		/*  first test overlap in the {x,y,z}-directions */

		/*  find min, max of the triangle each direction, and test for overlap in */

		/*  that direction -- this is equivalent to testing a minimal AABB around */

		/*  the triangle against the AABB */



		/* test in X-direction */
		FindMinMax(v0.x, v1.x, v2.x, min, max);
		//FINDMINMAX(v0[X], v1[X], v2[X], min, max);

		if (min > boxhalfsize.x || max < -boxhalfsize.x) return false;



		/* test in Y-direction */
		FindMinMax(v0.y, v1.y, v2.y, min, max);
		//FINDMINMAX(v0[Y], v1[Y], v2[Y], min, max);

		if (min > boxhalfsize.y || max < -boxhalfsize.y) return false;



		/* test in Z-direction */
		FindMinMax(v0.z, v1.z, v2.z, min, max);
		//FINDMINMAX(v0[Z], v1[Z], v2[Z], min, max);

		if (min > boxhalfsize.z || max < -boxhalfsize.z) return false;



		/* Bullet 2: */

		/*  test if the box intersects the plane of the triangle */

		/*  compute plane equation of triangle: normal*x+d=0 */
		normal = glm::cross(e0, e1);
		//CROSS(normal, e0, e1);

		// -NJMP- (line removed here)

		if (!planeBoxOverlap(normal, v0, boxhalfsize)) return false;	// -NJMP-


		return true;   /* box and triangle overlaps */
	}

}//namespace Moller

