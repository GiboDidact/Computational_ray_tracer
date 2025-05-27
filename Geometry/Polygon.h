#pragma once
#include "../pch.h"
#include "Segment.h"

class GraphicsManager;
class Polygon
{
public:
	Polygon();
	Polygon(std::vector<glm::vec2> data);

	void AddVertex(glm::vec2 pos);
	void RemoveVertex(int index);
	void RemoveVertexEar(int index);
	void splitEdge(int first_vertid);
	void joinEdge(int vertid);

	void Shift(glm::vec2 offset);
	void ShiftVertex(glm::vec2 offset, int idx);
	void Rotate(float theta);
	void Scale(float scale);
	glm::vec2 getCenter() const;

	void CalculateEars(bool shuffle = false, bool one_ear = false);
	glm::ivec4 FindEarInTriangulation();
	void CalculateMouth();
	float CalculateArea();
	float CalculateInternalAngles();
	bool CalculateCCW();
	glm::ivec2 CalculateDiagnol(Polygon P, bool only_ear = false, glm::ivec2 ear_exclusion = glm::ivec2(-1, -1));

	void DrawPolygon(GraphicsManager& g_manager, bool triangulation = false, bool highlight_vertex = false, glm::vec3 col = glm::vec3(-1, -1, -1),
		bool draw_edges = true, bool draw_vert = true) const;
	void HighlightEars(GraphicsManager& g_manager) const;
	void HighlightDiagnol(GraphicsManager& g_manager) const;
	void HighlightMouth(GraphicsManager& g_manager) const;

	bool IsConvex() const;
	bool isVertexReflex(int id) const;
	bool isValidPolygon() const;

	int hitVertex(glm::vec2 pos, float threshold) const;
	int hitBoundary(glm::vec2 pos, float threshold) const;
	int hitBoundaryNoVertex(glm::vec2 pos, float threshold) const;
	bool IsInside(glm::vec2 pos, glm::vec2 _dir = glm::vec2(-1, -1)) const;

	void TriangulateDiagnolSplitting(bool ears = false, bool shuffle = false);
	
	void FindDiagnol();
	void MakeCCW();
	int getreflexcount() const;

	void PrintVertexData() const;
	void PrintTriangleCount() const;
	void ExportPolygon(std::string file) const;

	int findVertexIndex(glm::vec2 v1) const;

	std::vector<glm::vec2>& V() { return vertices; }
	int gettrianglecount() const { return triangle_indices.size(); }
	int getvertexcount() const { return vertices.size(); }
	bool IsCCW() const { return CCW; }
	std::vector<glm::ivec3>& T() { return triangle_indices; }
	std::vector<glm::ivec3> getearindex() { return ears_indices; }
private:
	//a list of points, and their neighbors are by proxy
	std::vector<glm::vec2> vertices;

	//the indices of the triangulation of the polygon
	std::vector<glm::ivec3> triangle_indices;

	//the two ears of the polygon
	std::vector<glm::ivec3> ears_indices;
	
	//one diagnol calculated
	glm::ivec2 diagnol_indices;

	//one mouth calculated
	glm::ivec2 mouth_indices;

	//orientation CCW or CW
	bool CCW;
private:
	//debug variables
	bool print_debug = false;

	void TriangulateDiagnolSplittingRec(Polygon P, bool ears = false);


public:
	static Polygon GenerateRandomPolygon(int N)
	{
		bool found = false;
		Polygon p;
		int generated_times = 0;
		while (!found)
		{
			generated_times++;
			p = GenerateRandomPolygonInt(N);
			if (p.isValidPolygon())
				found = true;
		}
		p.MakeCCW();
		return p;
	}

	static Polygon GenerateRandomPolygonInt(int N)
	{
		//So we need to make sure its a valid simple polygon, its closed and it doesn't intersect itself
		Polygon p;

		//pick a random point
		glm::vec2 starting_point(0, 0);
		p.AddVertex(starting_point);
		int redo_amount = 0;
		for (int i = 1; i < N; i++)
		{
			bool found = false;
			//find a new edge
			while (!found)
			{
				redo_amount++;
				float length = Helper::GetRandomNumber(50, 100);
				glm::vec2 dir = glm::normalize(glm::vec2(Helper::GetRandomNumber(-1, 1), Helper::GetRandomNumber(-1, 1)));
				glm::vec2 new_pos = p.V()[i - 1] + dir * length;

				if (p.V().size() == 1)
				{
					p.AddVertex(new_pos);
					found = true;
				}
				else
				{
					bool intersects = false;
					for (int v = 0; v < p.V().size() - 1; v++)
					{
						int id1 = v;
						int id2 = v + 1;
						if (Segment::SegmentsIntersectNoBoundary(p.V()[id1], p.V()[id2], p.V()[i - 1], new_pos))
						{
							intersects = true;
							break;
						}
					}
					if (!intersects)
					{
						p.AddVertex(new_pos);
						found = true;
					}
				}
				//early outs
				if (N < 50 && redo_amount > 500)
				{
					p.vertices.clear();
					return p;
				}
				else if (N > 50 && redo_amount > 1000)
				{
					p.vertices.clear();
					return p;
				}
				else if (N > 250 && redo_amount > 2000)
				{
					p.vertices.clear();
					return p;
				}
			}
		}
		//std::cout << "redo amount: " << redo_amount << std::endl;

		return p;
	}

	static Polygon ImportPolygon(std::string filename)
	{
		//open the file and construct the polygon
		Polygon poly1;
		std::fstream polyfile;
		std::string name = "Game_Data/polygon_data/" + filename;
		polyfile.open(name, std::fstream::in);
		if (polyfile.is_open())
		{
			while (!polyfile.eof())
			{
				//read first number, comma, second number. nextline
				std::string x_data;
				std::getline(polyfile, x_data, ',');
				std::string y_data;
				std::getline(polyfile, y_data, '\n');
				//std::cout << x_data << " " << y_data << std::endl;
				float x_pos = std::stof(x_data);
				float y_pos = std::stof(y_data);

				poly1.AddVertex(glm::vec2(x_pos, y_pos));
			}

			if (polyfile.bad() || polyfile.fail())
			{
				std::cout << "polygon file opened but failed writing to file\n";
				return poly1;
			}
		}
		else
		{
			std::cout << "failed to open " << filename << std::endl;
		}

		return poly1;
	}
};