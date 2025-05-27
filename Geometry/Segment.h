#pragma once
#include "../pch.h"

namespace Segment
{
	//takes in 3 points, where first 1 is fixed
	static float DirectionSegments(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3)
	{
		//create two vectors at origin
		glm::vec2 vector1 = p2 - p1;
		glm::vec2 vector2 = p3 - p1;

		//now the cross product gives you a signed number +, -, or 0 if they are colinear
		return vector1.x*vector2.y - vector2.x*vector1.y;
	}

	static int TurnSegments(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3)
	{
		//positive value means it took a right
		//negative means it took a left

		float val = DirectionSegments(p1, p2, p3);
		if (val < 0)
			return 1;
		if (val > 0)
			return -1;
		return 0;
	}
	
	static bool onSegment(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3)
	{
		if ((std::min(p1.x, p2.x) <= p3.x && p3.x <= std::max(p1.x, p2.x)) && (std::min(p1.y, p2.y) <= p3.y  && p3.y <= std::max(p1.y, p2.y)))
			return true;

		return false;
	}

	static bool SegmentsIntersectNoBoundary(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3, glm::vec2 p4)
	{
		float d1 = DirectionSegments(p1, p2, p3);
		float d2 = DirectionSegments(p1, p2, p4);
		float d3 = DirectionSegments(p3, p4, p1);
		float d4 = DirectionSegments(p3, p4, p2);

		if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 < 0 && d4 > 0) || (d3 > 0 && d4 < 0)))
			return true;

		return false;
	}

	static bool SegmentsIntersect(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3, glm::vec2 p4)
	{
		float d1 = DirectionSegments(p1, p2, p3);
		float d2 = DirectionSegments(p1, p2, p4);
		float d3 = DirectionSegments(p3, p4, p1);
		float d4 = DirectionSegments(p3, p4, p2);

		if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 < 0 && d4 > 0) || (d3 > 0 && d4 < 0)))
			return true;
		else if (d1 == 0 && onSegment(p1, p2, p3))
			return true;
		else if (d2 == 0 && onSegment(p1, p2, p4))
			return true;
		else if (d3 == 0 && onSegment(p3, p4, p1))
			return true;
		else if (d4 == 0 && onSegment(p3, p4, p2))
			return true;

		return false;
	}

	/*https://iquilezles.org/articles/distfunctions2d/*/
	static float DistanceFromSegment(glm::vec2 p, glm::vec2 a, glm::vec2 b)
	{
		glm::vec2 pa = p - a;
		glm::vec2 ba = b - a;
		float h = std::clamp(glm::dot(pa, ba) / glm::dot(ba, ba), 0.0f, 1.0f);
		return glm::length(pa - ba * h);
	}

	/*https://iquilezles.org/articles/distfunctions2d/*/
	static float sdTriangle(glm::vec2 p, glm::vec2 p0, glm::vec2 p1, glm::vec2 p2)
	{
		glm::vec2 e0 = p1 - p0, e1 = p2 - p1, e2 = p0 - p2;
		glm::vec2 v0 = p - p0, v1 = p - p1, v2 = p - p2;
		glm::vec2 pq0 = v0 - e0 * std::clamp(glm::dot(v0, e0) / glm::dot(e0, e0), 0.0f, 1.0f);
		glm::vec2 pq1 = v1 - e1 * std::clamp(glm::dot(v1, e1) / glm::dot(e1, e1), 0.0f, 1.0f);
		glm::vec2 pq2 = v2 - e2 * std::clamp(glm::dot(v2, e2) / glm::dot(e2, e2), 0.0f, 1.0f);
		float s = glm::sign(e0.x * e2.y - e0.y * e2.x);
		glm::vec2 d = min(min(glm::vec2(dot(pq0, pq0), s * (v0.x * e0.y - v0.y * e0.x)),
			glm::vec2(dot(pq1, pq1), s * (v1.x * e1.y - v1.y * e1.x))),
			glm::vec2(dot(pq2, pq2), s * (v2.x * e2.y - v2.y * e2.x)));
		return -sqrt(d.x) * glm::sign(d.y);
	}

	static glm::vec2 SegmentToLine(glm::vec2 a, glm::vec2 b)
	{
		//make sure they are rightmost
		glm::vec2 p1 = a;
		glm::vec2 p2 = b;
		if (b.x < a.x) {
			p1 = b;
			p2 = a;
		}
		double m1 = (p2.y - p1.y) / (p2.x - p1.x);
		float t1 = -p1.x / (p2.x - p1.x); //t of ray at x=0
		double b1 = p1.y + t1 * (p2.y - p1.y);

		glm::vec2 answer(m1, b1);
		return answer;
	}

	//return: first tells you if it hit or not, next two give you x,y position if it
	static glm::vec3 LineToLineIntersectionPoint(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 d)
	{
		if (!SegmentsIntersect(a, b, c, d))
		{
			return glm::vec3(-1, 0, 0);
		}
		
		//calculate m1, m2, b1, b2
		glm::vec2 line1 = SegmentToLine(a, b);
		glm::vec2 line2 = SegmentToLine(c, d);
		std::cout << "m1: " << line1.x << " b1: " << line1.y << std::endl;
		std::cout << "m2: " << line2.x << " b2: " << line2.y << std::endl;
		
		//intersection point
		float x_pos = (line1.y - line2.y) / (line2.x - line1.x);
		float y_pos = line1.x * x_pos + line1.y;

		return glm::vec3(1, x_pos, y_pos);
	}

	static glm::vec3 LineToLineIntersectionPoint2(glm::vec2 a, glm::vec2 b, glm::vec2 c, glm::vec2 d)
	{
		if (!SegmentsIntersect(a, b, c, d))
		{
			return glm::vec3(-1, 0, 0);
		}

		double a1 = b.y - a.y;
		double b1 = a.x - b.x;
		double c1 = a1 * (a.x) + b1 * (a.y);

		double a2 = d.y - c.y;
		double b2 = c.x - d.x;
		double c2 = a2 * (c.x) + b2 * (c.y);

		double determinant = a1 * b2 - a2 * b1;

		if (determinant == 0)
		{
			//the lines are parallel
			return glm::vec3(-1, 999998, 0);
		}
		else
		{
			double x = (b2 * c1 - b1 * c2) / determinant;
			double y = (a1 * c2 - a2 * c1) / determinant;
			return glm::vec3(1,x, y);
		}
	}
}