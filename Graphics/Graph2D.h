#pragma once
#include "GraphicsManager.h"
class Graph2D
{
public:

	Graph2D(GraphicsManager& g_manager, int resolutionx, int resolutiony);
	~Graph2D();


	void DrawToTexture();
	void setbackgroundcolor(glm::vec3 col);
	void give_mouse_data(glm::vec2 mouse_position, int zoom, bool holding_left, bool first_left_click);
	void move_to_pos(glm::vec2 pos);
	void scale_cam(glm::vec2 scale);
	void setviewport(float x, float y, float width, float height);

	void plotpoint(int line, const std::pair<float, float>& point);
	//void plotpoints(const std::vector<std::pair<float, float>>& points);
	void plotpoints(int line, const std::vector<float>& points);
	void clearalllines();
	void clearline(int line);
	glm::vec2 mousetographpos(glm::vec2 pos);

	void createexponential();
	void createxsquared();
	void createsin(float A, float w);
	void createlinear();
	void createfunction(int line, float x_range, float x_precision, std::function<float(float)> graph_function);
	void createfunction(int line, float a, float b, float x_precision, std::function<float(float)> graph_function);
	
	void setlinecolor(glm::vec3 col, int line)
	{
		line_colors[line] = col;
	}

	void enableline(int line, bool enable)
	{
		active_lines[line] = enable;
	}

	GLuint gettexture() const { return color_attachment_blit_id; }
private:
	GraphicsManager& graphics_manager;

	glm::vec4 background_color;
	glm::vec4 viewport_data;
	glm::mat4 proj_m;
	glm::mat4 cam_m;

	glm::vec2 cam_center;
	glm::vec2 mouse_last_pos;

	GLsizei resolution_w;
	GLsizei resolution_h;
	GLfloat	cam_w;
	GLfloat cam_h;

	std::array<bool, 10> active_lines;
	std::array<glm::vec3, 10> line_colors;
	std::array<std::vector<float>,10> data_points; //(x,y) points paired in consecutive index. point 1 -->0,1, point 2-->2,3, point 3-->4,5
	std::vector<float> grid_points;
	std::vector<float> subgrid_points;
	std::vector<GraphicsManager::RenderText> rendertextqueue;

	GLuint fbo_id;
	GLuint color_attachment_id;
	std::array<GLuint,10> vao_points, vbo_points;
	GLuint vao_grids, vbo_grids;
	GLuint vao_subgrids, vbo_subgrids;
	GLuint fbo_blit_id;
	GLuint color_attachment_blit_id;

	void createdefaultgrids();
	void fillgriddata(float spacing, std::vector<float>& point_data, bool add_numbers, int exponent, float sub_space);
	void writegpubuffer(int line);
	void writegridbuffer();
	void writesubgridbuffer();
	uint32_t point_index(uint32_t point);
};

