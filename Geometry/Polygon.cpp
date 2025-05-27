#include "../pch.h"
#include "Polygon.h"
#include "../Graphics/GraphicsManager.h"

Polygon::Polygon()
{
	vertices.reserve(10);
	diagnol_indices = glm::ivec2(-1, -1);
	mouth_indices = glm::ivec2(-1, -1);
	CCW = false;
}
Polygon::Polygon(std::vector<glm::vec2> data)
{
	vertices.reserve(data.size());
	for (int i = 0; i < data.size(); i++)
	{
		vertices.push_back(data[i]);
	}

	diagnol_indices = glm::ivec2(-1, -1);
	mouth_indices = glm::ivec2(-1, -1);

	CCW = CalculateCCW();
}

void Polygon::AddVertex(glm::vec2 pos)
{
	vertices.push_back(pos);
	CCW = CalculateCCW();
}

void Polygon::RemoveVertex(int index)
{
	if (index <= vertices.size())
		vertices.erase(vertices.begin() + index);

	CCW = CalculateCCW();
}

void Polygon::RemoveVertexEar(int index)
{
	//remove triangle that was the ear, since its an eartip there should only be 1 triangle
	for (int t = 0; t < triangle_indices.size(); t++)
	{
		if (triangle_indices[t].x == index || triangle_indices[t].y == index || triangle_indices[t].z == index)
		{
			triangle_indices.erase(triangle_indices.begin() + t);
			break;
		}
	}

	if (index < vertices.size())
		vertices.erase(vertices.begin() + index);

	//now we have to update all the vertices in the triangle indices. Really if its before do nothing if its after subtract 1, shouldnt equal
	for (int i = 0; i < triangle_indices.size(); i++)
	{
		if (triangle_indices[i].x > index)
			triangle_indices[i].x -= 1;
		if (triangle_indices[i].y > index)
			triangle_indices[i].y -= 1;
		if (triangle_indices[i].z > index)
			triangle_indices[i].z -= 1;
	}

}

void Polygon::splitEdge(int first_vertid)
{
	//insert a new vertice between first_vertid and first_vertid + 1, and perturb a bit
	int idx2 = (first_vertid + 1) % vertices.size();
	glm::vec2 new_pos = vertices[first_vertid] * .5f + vertices[idx2] * .5f;
	glm::vec2 dir = glm::normalize(vertices[idx2] - vertices[first_vertid]);
	dir = glm::vec2(dir.y, -dir.x);
	new_pos += dir * glm::length(vertices[idx2] - vertices[first_vertid]) * .2f;

	vertices.insert(vertices.begin() + first_vertid + 1, new_pos);
}

void Polygon::joinEdge(int vertid)
{
	if (vertices.size() <= 3)
		return;
	vertices.erase(vertices.begin() + vertid);
}

void Polygon::Shift(glm::vec2 offset)
{
	for (int i = 0; i < vertices.size(); i++)
	{
		vertices[i] += offset;
	}
}

void Polygon::ShiftVertex(glm::vec2 offset, int idx)
{
	vertices[idx] += offset;
	CCW = CalculateCCW();
}

void Polygon::Rotate(float theta)
{
	//translate to origin, do rotation, translate back
	//the center will just be the average of all the  points
	glm::vec2 center = getCenter();

	for (int i = 0; i < vertices.size(); i++)
	{
		glm::vec2 pos = vertices[i];
		pos -= center;
		pos = glm::vec4(pos.x, pos.y, 0, 1) * glm::rotate(glm::mat4(1.0f), glm::radians(theta), glm::vec3(0, 0, 1));
		pos += center;

		vertices[i] = pos;
	}
}

void Polygon::Scale(float scale)
{
	//the center will just be the average of all the  points
	glm::vec2 center = getCenter();

	for (int i = 0; i < vertices.size(); i++)
	{
		glm::vec2 pos = vertices[i];
		pos -= center;
		pos = glm::vec4(pos.x, pos.y, 0, 1) * glm::scale(glm::mat4(1.0f), glm::vec3(scale, scale, 1));
		pos += center;

		vertices[i] = pos;
	}
}

glm::vec2 Polygon::getCenter() const
{
	glm::vec2 center(0, 0);
	for (int i = 0; i < vertices.size(); i++)
	{
		center += vertices[i];
	}
	center /= vertices.size();
	return center;
}

void Polygon::CalculateEars(bool shuffle, bool one_ear)
{
	//Corollary 1.9: Every polygon with more than three vertices has at least two ears

	ears_indices.clear();
	if (vertices.size() < 3)
		return;
	if (vertices.size() == 3)
	{
		ears_indices.push_back(glm::ivec3(0, 1, 2));
		return;
	}


	Polygon poly = *this;
	if (shuffle)
	{
		Helper::ShiftArray<glm::vec2>(poly.vertices.data(), poly.vertices.size(), Helper::GetRandomNumber(0, poly.vertices.size() - 1));
	}

	//Ear one
	glm::ivec2 ids = CalculateDiagnol(poly, true);
	glm::ivec2 ids_save = ids;
	ids.x = findVertexIndex(poly.V()[ids.x]);
	ids.y = findVertexIndex(poly.V()[ids.y]);
	int id_middle = (ids.x + 1) % vertices.size();
	ears_indices.push_back(glm::ivec3(ids.x, id_middle, ids.y));

	if (one_ear)
		return;

	//Ear two
	glm::ivec2 ids2 = CalculateDiagnol(poly, true, ids_save);
	ids2.x = findVertexIndex(poly.V()[ids2.x]);
	ids2.y = findVertexIndex(poly.V()[ids2.y]);
	id_middle = (ids2.x + 1) % vertices.size();
	ears_indices.push_back(glm::ivec3(ids2.x, id_middle, ids2.y));
}

glm::ivec4 Polygon::FindEarInTriangulation()
{
	if (triangle_indices.empty())
		return glm::ivec4(-1, -1, -1, -1);

	//just go through all the triangles and find one that is consecutive
	glm::ivec4 ear(-1, -1, -1, -1);
	for (int t = 0; t < triangle_indices.size(); t++)
	{
			int id1 = triangle_indices[t].x;
			int id2 = triangle_indices[t].y;
			int id3 = triangle_indices[t].z;


			//from id1 it either has to be double front, double back, or 1 from each
			//search for front
			int front = -1;
			if ((id1 + 1) % vertices.size() == id2)
				front = id2;
			else if ((id1 + 1) % vertices.size() == id3)
				front = id3;
			//if we find front see if other one is in front of that
			if (front == id2 && (id2 + 1) % vertices.size() == id3)
			{
				ear = glm::ivec4(id1, id2, id3, t);
				break;
			}
			if (front == id3 && (id3 + 1) % vertices.size() == id2)
			{
				ear = glm::ivec4(id1, id3, id2, t);
				break;
			}
			//see if the back
			if (front == id2 && (id1 + vertices.size() - 1) % vertices.size() == id3)
			{
				ear = glm::ivec4(id3, id1, id2, t);
				break;
			}
			if (front == id3 && (id1 + vertices.size() - 1) % vertices.size() == id2)
			{
				ear = glm::ivec4(id2, id1, id3, t);
				break;
			}

			int back = -1;
			if ((id1 + vertices.size() - 1) % vertices.size() == id2)
				back = id2;
			else if ((id1 + vertices.size() - 1) % vertices.size() == id3)
				back = id3;
			if (back == id2 && (id2 + vertices.size() - 1) % vertices.size() == id3)
			{
				ear = glm::ivec4(id3, id2, id1, t);
				break;
			}
			if (back == id3 && (id3 + vertices.size() - 1) % vertices.size() == id2)
			{
				ear = glm::ivec4(id2, id3, id1, t);
				break;
			}
	}

	return ear;
}

void Polygon::CalculateMouth()
{
	//One-Mouth Theorem: Except for convex polygons, every simple polygon has at least one mouth
	//If a reflex vertex is not a mouth, then it contains a convex vertex
	//converse: if the triangle of a relfex vertex does not contain a convex vertex, then it is a mouth
	//pretty much it either containts a convex with anything or if it doesnt have a convex it cant have a reflex or anything else
	
	mouth_indices = glm::ivec2(-1, -1);
	if (IsConvex())
		return;

	if (vertices.size() <= 3)
		return;

	Polygon poly = *this;
	Helper::ShiftArray<glm::vec2>(poly.vertices.data(), poly.vertices.size(), Helper::GetRandomNumber(0, poly.vertices.size() - 1));
	bool found = false;
	for (int i = 0; i < poly.V().size(); i++)
	{
		//find reflex vertex
		if (poly.isVertexReflex(i))
		{
			//see if anything is inside the triangle
			bool intersection = false;
			bool one_concave = false; //sanity check, there has to always be at least 1 concave out of all intersections

			int id_behind = (i + poly.V().size() - 1) % poly.V().size();
			int id_forward = (i + 1) % poly.V().size();

			Polygon tri1;
			tri1.AddVertex(poly.V()[id_behind]);
			tri1.AddVertex(poly.V()[i]);
			tri1.AddVertex(poly.V()[id_forward]);
			//technically we have to check that it doesnt hit a convex vertex, but I don't think it matters if there was a vertex hit at least 1 has to be convex
			for (int j = 0; j < poly.V().size() - 3; j++)
			{
				int front_id = (id_forward + 1 + j) % poly.V().size();
				if (tri1.IsInside(poly.V()[front_id]))
				{
					if (!poly.isVertexReflex(front_id));
					one_concave = true;
					intersection = true;
				}
			}
			if (!intersection)
			{
				mouth_indices = glm::ivec2(id_behind, id_forward);
				found = true;
				break;
			}
			else
			{
				if (one_concave == false)
				{
					std::cout << "all the points inside the mouth are reflex??\n";
				}
			}
		}
	}
	if (found)
	{
		mouth_indices.x = findVertexIndex(poly.V()[mouth_indices.x]);
		mouth_indices.y = findVertexIndex(poly.V()[mouth_indices.y]);
	}
	else
	{
		std::cout << "couldn't find mouth????\n";
	}
}

float Polygon::CalculateArea()
{
	//The area can be derived from Green's Theorem

	float sum = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		int id1 = i;
		int id2 = (i + vertices.size() - 1) % vertices.size();
		sum += (vertices[id1].x * vertices[id2].y - vertices[id2].x * vertices[id1].y);
	}
	sum = abs(sum);
	sum /= 2.0f;

	return sum;
}

float Polygon::CalculateInternalAngles()
{
	//sum of interior angles is always pi*(n-2)

	float total_angle = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		int id1 = i;
		int id2 = (i + 1) % vertices.size();
		int id3 = (i + 2) % vertices.size();

		glm::vec2 v1 = glm::normalize(vertices[id1] - vertices[id2]);
		glm::vec2 v2 = glm::normalize(vertices[id3] - vertices[id2]);

		if (isVertexReflex(id2))
			total_angle += (2 * 3.14) - acos(glm::dot(v1, v2));
		else
			total_angle += acos(glm::dot(v1, v2));
	}

	return total_angle;
}

//we know how to find if the boundary takes a left or right turn, and convex/reflex vertices have different directions, CCW or CW has them vice versa'd
//in  CCW: convex-->left turn, reflex-->right turn   in CW: convex-->right turn, reflex-->left turn
//find the bottommost/rightmost point, you know this has to be a convex. Now if its right its CW, if its left its CCW
bool Polygon::CalculateCCW()
{
	if (vertices.size() <= 2)
		return false;

	int THE_ID = -1;
	std::vector<int> bottom_most_ids;
	int min_value = 1000000.0f;
	//find min value
	for (int i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].y < min_value)
			min_value = vertices[i].y;
	}
	//get all points within min value
	for (int i = 0; i < vertices.size(); i++)
	{
		if (abs(vertices[i].y - min_value) <= 1)
			bottom_most_ids.push_back(i);
	}

	if (bottom_most_ids.size() == 0)
	{
		std::cout << "no bottommostid??\n";
		return false;
	}
	else if (bottom_most_ids.size() == 1)
	{
		THE_ID = bottom_most_ids[0];
	}
	else
	{
		//shared points, find the rightmost one out of them
		int right_value = -1000000.0;
		for (int i = 0; i < bottom_most_ids.size(); i++)
		{
			if (vertices[bottom_most_ids[i]].x >= right_value)
			{
				THE_ID = bottom_most_ids[i];
			}
		}
	}
	if (THE_ID == -1) {
		std::cout << "CCW error, no id found\n";
		return false;
	}

	int id_before = (THE_ID + vertices.size() - 1) % vertices.size();
	int id_after = (THE_ID + 1) % vertices.size();

	if (Segment::TurnSegments(vertices[id_before], vertices[THE_ID], vertices[id_after]) == -1)
	{
		return true;
	}
	return false;
}

glm::ivec2 Polygon::CalculateDiagnol(Polygon P, bool only_ear, glm::ivec2 ear_exclusion)
{
	if (P.V().size() <= 3)
	{
		return glm::vec2(-1, -1);
	}

	//go through every nonadjacent pairs(all possible diagnols) and see if it intersects any of the edges of the polygon(except at the endpoint)
	bool diagnol_found = false;
	glm::ivec2 diagnol_vertices(0, 0);
	for (int i = 0; i < P.V().size(); i++)
	{
		for (int j = i + 2; j < P.V().size(); j++)
		{
			glm::vec2 v1 = P.V()[i];
			glm::vec2 v2 = P.V()[j];
			//neighbors, ignore
			if (i == 0 && j == P.V().size() - 1)
				break;
			if (only_ear)
			{
				//we only care about diagnols that are 1 away from eachother
				if (j - i != 2)
					break;

				//if we want to find an ear thats not this one (every polygon has at least 2 ears)
				//**we don't want ears that are intersecting(right next to eachother)
				if (ear_exclusion.x != -1 && ear_exclusion.y != -1)
				{
					if (i == ear_exclusion.x && j == ear_exclusion.y)
						break;

					int idx_behind = (ear_exclusion.x + P.V().size() - 1) % P.V().size();
					int idx_infront = (ear_exclusion.x + 1) % P.V().size();
					int idy_infront = (ear_exclusion.y + 1) % P.V().size();
					if (i == idx_behind && j == idx_infront)
						break;
					if (i == idx_infront && j == idy_infront)
						break;
				}
			}

			//see if this diagnol intersects with the boundary
			bool intersects = false;
			for (int x = 0; x < P.V().size(); x++)
			{
				int idx = x + 1;
				if (x == P.V().size() - 1)
					idx = 0;
				glm::vec2 e1 = P.V()[x];
				glm::vec2 e2 = P.V()[idx];
				if (Segment::SegmentsIntersectNoBoundary(v1, v2, e1, e2))
				{
					intersects = true;
					break;
				}
			}
			if (!intersects)
			{
				//so we know it doesnt intersect any of the edges, but it can still be completely outside or completely inside
				//we pick 1 point maybe 2 to be safe and see if its inside or out. Since it doesnt cross the boundary all points are on the same side
				//by jordan curve theorem

				//pick a number between 0 and 1 to chooose a point along the line At + B(1-t)
				float t1 = .2f;
				float t2 = .7f;
				glm::vec2 p1 = P.V()[i] * t1 + P.V()[j] * (1.0f - t1);
				glm::vec2 p2 = P.V()[i] * t2 + P.V()[j] * (1.0f - t2);

				bool inside1 = P.IsInside(p1);
				bool inside2 = P.IsInside(p2);

				if (inside1 && inside2)
				{
					diagnol_found = true;
					diagnol_vertices = glm::ivec2(i, j);
					break;
				}

			}
		}
		if (diagnol_found)
			break;
	}

	if (!diagnol_found)
		std::cout << "didn't find a diagnol???\n";

	return diagnol_vertices;
}

void Polygon::DrawPolygon(GraphicsManager& g_manager, bool triangulation, bool highlight_vertex, glm::vec3 col, bool draw_edges, bool draw_vert) const
{
	glm::vec3 interior_color(50, 100, 100);
	if (col != glm::vec3(-1, -1, -1))
		interior_color = col;
	glm::vec3 vertex_color(0, 0, 0);
	glm::vec3 edge_color(0, 0, 0);

	//INTERIOR
	if (vertices.size() >= 3)
	{
		for (int i = 0; i < triangle_indices.size(); i++)
		{
			int idx1 = triangle_indices[i].x;
			int idx2 = triangle_indices[i].y;
			int idx3 = triangle_indices[i].z;
			glm::vec2 first_pos = glm::vec2(vertices[idx1].x, vertices[idx1].y);
			glm::vec2 second_pos = glm::vec2(vertices[idx2].x, vertices[idx2].y);
			glm::vec2 third_pos = glm::vec2(vertices[idx3].x, vertices[idx3].y);

			std::array<glm::vec2, 3> points = { first_pos, second_pos, third_pos };
			g_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::DYN_TRIANGLE, SHADER_TYPE::MAIN, glm::mat4(1.0f), interior_color, glm::vec4(0,0,0,0), points });
		}
	}

	//EDGES
	if (triangulation)
	{
		for (int i = 0; i < triangle_indices.size(); i++)
		{
			glm::vec2 v1 = glm::vec2(vertices[triangle_indices[i].x].x, vertices[triangle_indices[i].x].y);
			glm::vec2 v2 = glm::vec2(vertices[triangle_indices[i].y].x, vertices[triangle_indices[i].y].y);
			glm::vec2 v3 = glm::vec2(vertices[triangle_indices[i].z].x, vertices[triangle_indices[i].z].y);

			g_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), edge_color, glm::vec4(v1.x,v1.y,v2.x,v2.y) });
			g_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), edge_color, glm::vec4(v2.x,v2.y,v3.x,v3.y) });
			g_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), edge_color, glm::vec4(v1.x,v1.y,v3.x,v3.y) });
		}
	}
	else
	{
		if (draw_edges)
		{
			for (int i = 0; i < vertices.size(); i++)
			{
				int indx = i + 1;
				if (i == vertices.size() - 1)
					indx = 0;

				glm::vec2 first_pos = glm::vec2(vertices[i].x, vertices[i].y);
				glm::vec2 second_pos = glm::vec2(vertices[indx].x, vertices[indx].y);

				glm::vec4 positions(first_pos.x, first_pos.y, second_pos.x, second_pos.y);

				g_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), edge_color, positions });
			}
		}
	}

	//VERTICES
	if (draw_vert)
	{
		glm::mat4 mat(1.0f);
		for (int i = 0; i < vertices.size(); i++)
		{
			glm::vec2 first_pos = glm::vec2(vertices[i].x, vertices[i].y);

			mat = glm::translate(glm::mat4(1.0f), glm::vec3(first_pos.x, first_pos.y, 0));
			mat = glm::scale(mat, glm::vec3(5, 5, 1));
			glm::vec3 colz = vertex_color;
			if (highlight_vertex && i == 0)
				colz = glm::vec3(255, 0, 0);
			if (highlight_vertex && i == 1)
				colz = glm::vec3(0, 255, 0);
			if (highlight_vertex && i == 2)
				colz = glm::vec3(0, 0, 255);

			g_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, colz });
		}
	}
}

void Polygon::HighlightEars(GraphicsManager& g_manager) const
{
	for (int i = 0; i < ears_indices.size(); i++)
	{
		glm::vec2 first_pos = glm::vec2(vertices[ears_indices[i].x].x, vertices[ears_indices[i].x].y);
		glm::vec2 second_pos = glm::vec2(vertices[ears_indices[i].y].x, vertices[ears_indices[i].y].y);
		glm::vec2 third_pos = glm::vec2(vertices[ears_indices[i].z].x, vertices[ears_indices[i].z].y);

		std::array<glm::vec2, 3> points = { first_pos, second_pos, third_pos };
		g_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::DYN_TRIANGLE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(255,255,0), glm::vec4(0,0,0,0), points });
	}
}

void Polygon::HighlightDiagnol(GraphicsManager& g_manager) const
{
	if (diagnol_indices.x != -1 && diagnol_indices.y != -1)
	{
		glm::vec2 v1 = glm::vec2(vertices[diagnol_indices.x].x, vertices[diagnol_indices.x].y);
		glm::vec2 v2 = glm::vec2(vertices[diagnol_indices.y].x, vertices[diagnol_indices.y].y);

		g_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(255,0,100), glm::vec4(v1.x,v1.y,v2.x,v2.y) });
	}
}

void Polygon::HighlightMouth(GraphicsManager& g_manager) const
{
	if (mouth_indices.x != -1 && mouth_indices.y != -1)
	{
		glm::vec2 v1 = glm::vec2(vertices[mouth_indices.x].x, vertices[mouth_indices.x].y);
		glm::vec2 v2 = glm::vec2(vertices[mouth_indices.y].x, vertices[mouth_indices.y].y);

		g_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(0,255,0), glm::vec4(v1.x,v1.y,v2.x,v2.y) });
	}
}


bool Polygon::IsConvex() const
{
	//if we travel along each edge and we only take left or only take right turns then its convex
	int first_val = 0;
	bool chosen = false;
	for (int i = 0; i < vertices.size(); i++)
	{
		int id1 = i;
		int id2 = (i + 1) % vertices.size();
		int id3 = (i + 2) % vertices.size();

		int val = Segment::TurnSegments(vertices[id1], vertices[id2], vertices[id3]);
		if (!chosen)
		{
			if (val != 0)
			{
				first_val = val;
				chosen = true;
			}
		}
		else
		{
			//if its colinear just keep going
			if (val == 0)
				continue;

			if ((val > 0 && first_val < 0) || (val < 0 && first_val > 0))
				return false;
		}
	}
	return true;
}

bool Polygon::isVertexReflex(int id) const
{
	int id_behind = (id + vertices.size() - 1) % vertices.size();
	int id_forward = (id + 1) % vertices.size();

	//works with CCW and CW
	int value = 0;
	if (IsCCW())
		value = 1;
	else
		value = -1;

	if (Segment::TurnSegments(vertices[id_behind], vertices[id], vertices[id_forward]) == value)
	{
		return true;
	}
	return false;
}

bool Polygon::isValidPolygon() const
{
	//check whether its closed, and whether it has self-intersection
	if (vertices.size() <= 2)
		return false;

	//go through every edge and see if theres an intersection with any other edge on the boundary
	for (int i = 0; i < vertices.size() - 1; i++)
	{
		int e1_id1 = i;
		int e1_id2 = i + 1;

		for (int j = i + 1; j < vertices.size(); j++)
		{
			int e2_id1 = j;
			int e2_id2 = j + 1;
			if (j == vertices.size() - 1)
				e2_id2 = 0;

			if (Segment::SegmentsIntersectNoBoundary(vertices[e1_id1], vertices[e1_id2], vertices[e2_id1], vertices[e2_id2]))
			{
				return false;
			}
		}
	}
	return true;
}

int Polygon::hitVertex(glm::vec2 pos, float threshold) const
{
	for (int i = 0; i < vertices.size(); i++)
	{
		if (glm::distance(pos, vertices[i]) <= threshold)
		{
			return i;
		}
	}
	return -1;
}

int Polygon::hitBoundary(glm::vec2 pos, float threshold) const
{
	for (int i = 0; i < vertices.size(); i++)
	{
		int idx2 = i + 1;
		if (i == vertices.size() - 1) {
			idx2 = 0;
		}
		if (Segment::DistanceFromSegment(pos, vertices[i], vertices[idx2]) <= threshold)
		{
			return i;
		}
	}
	return -1;
}

int Polygon::hitBoundaryNoVertex(glm::vec2 pos, float threshold) const
{
	if (hitVertex(pos, threshold) != -1)
		return -1;

	int id = hitBoundary(pos, threshold);
	if (id != -1)
		return id;

	return -1;
}

bool Polygon::IsInside(glm::vec2 pos, glm::vec2 _dir) const
{
	if (vertices.size() <= 2)
	{
		return false;
	}
	//Jordan Curve Theorem says we can create a ray and count the intersections with the edges which partitions it into odd or even which correspond
	//to whether its inside or outside of the polygon. (even=outisde, odd=inside)

	//first see if its on the boundary, if its near the boundary we'll just force it to be inside. This helps against degenerate cases
	float eps2 = 0.0001f;
	if (hitBoundary(pos, eps2) != -1)
		return true;

	//choose a fixed direction that isn't the same as any edge, edges are finite so we can do this
	bool found_direction = false;
	glm::vec2 dir;
	float eps3 = .00001f;
	if (_dir.x == -1 && _dir.y == -1)
	{
		while (!found_direction)
		{
			dir = glm::vec2(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1));
			dir = glm::normalize(dir);

			bool too_close = false;
			for (int i = 0; i < vertices.size() - 1; i++)
			{
				glm::vec2 edgedir = glm::normalize(vertices[i + 1] - vertices[i]);
				if (1 - abs(glm::dot(edgedir, dir)) <= eps3)
				{
					too_close = true;
					break;
				}
			}
			if (too_close == false)
				found_direction = true;
		}
	}
	else
	{
		dir = _dir;
	}

	//create a ray with the input pos and direction, and count the number of intersections
	glm::vec2 pos2 = pos + dir * 100000.0f;

	int intersection_count = 0;
	int vertexhit = 0;
	std::vector<glm::ivec2> hitedges;
	for (int i = 0; i < vertices.size(); i++)
	{
		int idx2 = i + 1;
		if (i == vertices.size() - 1) {
			idx2 = 0;
		}
		if (Segment::SegmentsIntersect(pos, pos2, vertices[i], vertices[idx2]))
		{
			hitedges.push_back(glm::ivec2(i, idx2));
			intersection_count++;
		}
	}

	//deal with cases where we are on a vertex
	float eps = 0.000001f;
	for (int i = 0; i < hitedges.size(); i++)
	{
		//see if its close to a vertex
		if (abs(Segment::DirectionSegments(pos, pos2, vertices[hitedges[i].y])) <= eps)
		{
			int wrapid = (hitedges[i].y + 1) % vertices.size();
			int wrapid2 = (hitedges[i].y + vertices.size() - 1) % vertices.size();
			//if it is see if the other edge was hit
			if (Segment::SegmentsIntersect(pos, pos2, vertices[hitedges[i].y], vertices[wrapid]))
			{
				//see if they lie on different sides
				float dir1 = Segment::DirectionSegments(pos, pos2, vertices[wrapid2]);
				float dir2 = Segment::DirectionSegments(pos, pos2, vertices[wrapid]);
				if ((dir1 > 0 && dir2 < 0) || (dir1 < 0 && dir2 > 0))
				{
					//they got counted twice, get rid of one
					vertexhit++;
				}
				else if (dir1 == 0 && dir2 == 0)
				{
					vertexhit++;
				}
			}
		}
	}

	intersection_count -= vertexhit;

	if (intersection_count % 2 == 1)
		return true;

	return false;
}

void Polygon::FindDiagnol()
{
	if (vertices.size() <= 3)
		return;

	Polygon poly = *this;
	Helper::ShiftArray<glm::vec2>(poly.vertices.data(), poly.vertices.size(), Helper::GetRandomNumber(0, poly.vertices.size() - 1));
	diagnol_indices = CalculateDiagnol(poly);

	diagnol_indices.x = findVertexIndex(poly.V()[diagnol_indices.x]);
	diagnol_indices.y = findVertexIndex(poly.V()[diagnol_indices.y]);
}

void Polygon::MakeCCW()
{
	if (!IsCCW())
	{
		std::reverse(vertices.begin(), vertices.end());
		TriangulateDiagnolSplitting(true);
		CCW = CalculateCCW();
	}
}

int Polygon::getreflexcount() const
{
	int count = 0;
	for (int i = 0; i < vertices.size();i++)
	{
		if (isVertexReflex(i))
			count++;
	}
	return count;
}

void Polygon::PrintVertexData() const
{
	std::cout << "----vertices----\n";
	for (int i = 0; i < vertices.size(); i++)
	{
		std::cout << "vertex " << i << " (" << vertices[i].x << ", " << vertices[i].y << ")\n";
	}
}
void Polygon::PrintTriangleCount() const
{
	//Theorem 1.8: Every triangulation of a polygon P with n vertices has n-2 triangles and n-3 diagnols
	std::cout << "Vertices: " << vertices.size() << ", triangles: " << triangle_indices.size() << std::endl;
}

void Polygon::ExportPolygon(std::string file) const
{
	//open a file in Game_Data/polygon_data and just write 2 number seperated by comma for vertex on line
	std::fstream polyfile;
	std::string name = "Game_Data/polygon_data/" + file;
	polyfile.open(name, std::fstream::out);
	if (polyfile.is_open())
	{
		for (int i = 0; i < vertices.size(); i++)
		{
			std::string vertex_string = std::to_string(vertices[i].x) + "," + std::to_string(vertices[i].y);
			if (i != vertices.size() - 1)
				vertex_string += "\n";
			polyfile.write(vertex_string.data(), vertex_string.size());
		}

		if (polyfile.bad() || polyfile.fail())
		{
			std::cout << "polygon file opened but failed writing to file\n";
			return;
		}
	}
	else
	{
		std::cout << "failed to open " << name << std::endl;
	}
}

int Polygon::findVertexIndex(glm::vec2 v1) const
{
	for (int i = 0; i < vertices.size(); i++)
	{
		if (v1 == vertices[i])
			return i;
	}

	std::cout << "couldnt find vertex index error\n";
	return -1;
}

void Polygon::TriangulateDiagnolSplitting(bool ears, bool shuffle)
{
	//lemma 1.3: Every polygon with more than three vertices has a diagnol
	//Theorem 1.4: Every polygon has a triangulation
	//Corollary 1.9: Every polygon with more than three vertices has at least two ears
	if (print_debug)
		std::cout << "----TriangulateDiagnolSplitting----\n";
	if (vertices.size() <= 2)
	{
		return;
	}

	triangle_indices.clear();
	if (shuffle)
	{
		Polygon poly = *this;
		//poly.PrintVertexData();
		Helper::ShiftArray<glm::vec2>(poly.vertices.data(), poly.vertices.size(), Helper::GetRandomNumber(0, poly.vertices.size() - 1));
		//poly.PrintVertexData();
		TriangulateDiagnolSplittingRec(poly, ears);
	}
	else
	{
		TriangulateDiagnolSplittingRec(*this, ears);
	}

	if (triangle_indices.size() != vertices.size() - 2)
	{
		if (print_debug)
			std::cout << "wrong triangle count after triangulation!\n";
	}
}

void Polygon::TriangulateDiagnolSplittingRec(Polygon P, bool ears)
{
	if (print_debug)
		std::cout << "split call\n";
	if (P.V().size() <= 3)
	{
		if (P.V().size() != 3) {
			std::cout << "error not 3 vertices\n";
			return;
		}
		int id1 = findVertexIndex(P.V()[0]);
		int id2 = findVertexIndex(P.V()[1]);
		int id3 = findVertexIndex(P.V()[2]);
		triangle_indices.push_back(glm::ivec3(id1, id2, id3));
		if (print_debug)
			std::cout << id1 << " " << id2 << " " << id3 << " ending\n";
		return;
	}

	glm::ivec2 indices = CalculateDiagnol(P, ears);
	if (indices.x == 0 && indices.y == 0) {
		if (print_debug)
			std::cout << "couldn't find diagnol\n";
		return;
	}

	if (print_debug)
	{
		std::cout << "diagnol: (" << indices.x << ", " << indices.y << ")\n";
	}
	//start at id1 then go right until you hit id2.
	Polygon poly1;
	for (int i = 0; i < P.V().size(); i++)
	{
		int wrapid = (indices.x + i) % P.V().size();

		poly1.AddVertex(P.V()[wrapid]);
		if (print_debug)
			std::cout << "wrapid1: " << findVertexIndex(P.V()[wrapid]) << ", ";
		if (wrapid == indices.y)
			break;
	}
	if (print_debug)
		std::cout << std::endl;
	//start at id1 then go left until you hit id2
	Polygon poly2;
	for (int i = 0; i < P.V().size(); i++)
	{
		int wrapid = (indices.x + P.V().size() - i) % P.V().size();

		poly2.AddVertex(P.V()[wrapid]);
		if (print_debug)
			std::cout << "wrapid2: " << findVertexIndex(P.V()[wrapid]) << ", ";
		if (wrapid == indices.y)
			break;
	}
	if (print_debug)
		std::cout << std::endl;


	TriangulateDiagnolSplittingRec(poly1);
	TriangulateDiagnolSplittingRec(poly2);
}