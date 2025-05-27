#include "../pch.h"
#include "Graph2D.h"
#include <format>
#include <sstream>
#include "../Util/HelperFunctions.h"

Graph2D::Graph2D(GraphicsManager& g_manager, int resolutionx, int resolutiony) : graphics_manager(g_manager), background_color(0, 0, 0, 1), viewport_data(0, 0, 100, 100),
  cam_center(0, 0), proj_m(1.0f), cam_m(1.0f)
{
	active_lines.fill(true);
	line_colors.fill(glm::vec3(255, 0, 0));

	//create FBO, color_attachment
	resolution_w = resolutionx;
	resolution_h = resolutiony;
	cam_w = 20;
	cam_h = 20;

	//resolve fbo
	glGenTextures(1, &color_attachment_blit_id);
	glGenFramebuffers(1, &fbo_blit_id);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo_blit_id);
	glBindTexture(GL_TEXTURE_2D, color_attachment_blit_id);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, resolution_w, resolution_h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color_attachment_blit_id, 0);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	
	//multisampled fbo
	glGenTextures(1, &color_attachment_id);
	glGenFramebuffers(1, &fbo_id);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);

	glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, color_attachment_id);
	glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 8, GL_RGB, resolution_w, resolution_h, GL_TRUE);
	//glTexParameteri(GL_TEXTURE_2D_MULTISAMPLE, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	//glTexParameteri(GL_TEXTURE_2D_MULTISAMPLE, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, color_attachment_id, 0);

	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
	{
		std::cout << "Frame buffer failed to be complete on Graph2D!\n";
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	graphics_manager.checkError("framebuffer graph2d");

	proj_m = glm::ortho(-(cam_w / 2.0f), cam_w / 2.0f, -(cam_h / 2.0f), cam_h / 2.0f, -1.0f, 1.0f);
	cam_m = glm::mat4(1.0f);

	//create vao/vbo gpu buffer
	for (int i = 0; i < vao_points.size(); i++)
	{
		glGenVertexArrays(1, &vao_points[i]);
		glBindVertexArray(vao_points[i]);
		glGenBuffers(1, &vbo_points[i]);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_points[i]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
		glEnableVertexAttribArray(0);
	}

	//create grid gpu buffer
	glGenVertexArrays(1, &vao_grids);
	glBindVertexArray(vao_grids);
	glGenBuffers(1, &vbo_grids);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_grids);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
	glEnableVertexAttribArray(0);

	//create subgrid gpu buffer
	glGenVertexArrays(1, &vao_subgrids);
	glBindVertexArray(vao_subgrids);
	glGenBuffers(1, &vbo_subgrids);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_subgrids);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
	glEnableVertexAttribArray(0);

	createdefaultgrids();

	//std::cout <<"cam dimensions: "<< cam_w << " " << cam_h << std::endl;
}

Graph2D::~Graph2D()
{
	glDeleteTextures(1, &color_attachment_id);
	glDeleteFramebuffers(1, &fbo_id);

	glDeleteTextures(1, &color_attachment_blit_id);
	glDeleteFramebuffers(1, &fbo_blit_id);
	for (int i = 0; i < vao_points.size(); i++)
	{
		glDeleteBuffers(1, &vao_points[i]);
		glDeleteBuffers(1, &vbo_points[i]);
	}
	glDeleteBuffers(1, &vao_grids);
	glDeleteBuffers(1, &vbo_grids);

	glDeleteBuffers(1, &vao_subgrids);
	glDeleteBuffers(1, &vbo_subgrids);
}

void Graph2D::setbackgroundcolor(glm::vec3 col)
{
	background_color = glm::vec4(col.x,col.y,col.z,1);
}

glm::vec2 Graph2D::mousetographpos(glm::vec2 pos)
{
	glm::vec2 mouse_world_pos = graphics_manager.WindowToWorldGeneral(pos, graphics_manager.getWidth(), graphics_manager.getHeight(),
		viewport_data, cam_w, cam_h, cam_center.x, cam_center.y);

	return mouse_world_pos;
}

void Graph2D::move_to_pos(glm::vec2 pos)
{
	cam_center = pos;
	createdefaultgrids();
}

void Graph2D::scale_cam(glm::vec2 scale)
{
	cam_w = scale.x;
	cam_h = scale.y;
	createdefaultgrids();
}

//call every frame
void Graph2D::give_mouse_data(glm::vec2 mouse_position, int zoom, bool holding_left, bool first_left_click)
{
	//the camera position is whertever it is, but if we click and drag our mouse we move the camera the same world distance as our cursor moves
	//so if your cursor starts at (200,300) window position, and now we drag it to (150,300), well we need to move the camera -50 x direction in world
	//really your mouse is in some spot on the camera, so percentage-->cam-->world

	//graph2D can technically be anywhere in the viewport, so we have to account for that, if its covering the whole window then its just window
	//but if its like in some smaller viewport of the window you have to take that into account
	
	glm::vec2 mouse_world_pos = graphics_manager.WindowToWorldGeneral(mouse_position, graphics_manager.getWidth(), graphics_manager.getHeight(),
		viewport_data, cam_w, cam_h, cam_center.x, cam_center.y);

	//if mouse is out of camera dont do anything
	if (!Helper::pointAABBCollision(glm::vec4(cam_center.x, cam_center.y, cam_w, cam_h), mouse_world_pos))
	{
		return;
	}

	//std::cout << mouse_world_pos.x << " " << mouse_world_pos.y << std::endl;
	if (holding_left)
	{
		glm::vec2 mouse_world_delta = mouse_world_pos - mouse_last_pos;
		cam_center -= mouse_world_delta;

		createdefaultgrids();
	}
	
	if(first_left_click)
		mouse_last_pos = mouse_world_pos;
	
	//zoom is a perctange of the current cam_width and cam_height
	float percentage = .05; //5%
	if (zoom == 1)
	{
		cam_w = cam_w * (1 - percentage);
		cam_h = cam_h * (1 - percentage);

		cam_center += (mouse_world_pos - cam_center) * percentage;
	}
	else if (zoom == -1)
	{
		cam_w = cam_w * (1 + percentage);
		cam_h = cam_h * (1 + percentage);

		cam_center -= (mouse_world_pos - cam_center)* percentage;
	}
	
	//update grids if we zoom in
	if (zoom == 1 || zoom == -1)
	{
		createdefaultgrids();
	}
}

void Graph2D::createdefaultgrids()
{
	//1/2/5 * 10^x   -->1,2,5,10,20,50,100,200,500,etc.
	//1-->(0,10), 2-->(10,20), 5-->(20,50), 10-->(50,100),20(200),50(500),100(1000)
	int theexponent = std::floorf(log10(cam_w));
	int thebase = 1;
	float spacing = 0.0f;
	if (cam_w >= 1)
	{
		if ((log10(cam_w) - theexponent) <= .3 )
		{
			spacing = 1 * std::powf(10, theexponent - 1);
			thebase = 1;
		}
		else if ((log10(cam_w) - theexponent) <= .7)
		{
			spacing = 2 * std::powf(10, theexponent - 1);
			thebase = 2;
		}
		else
		{
			spacing = 5 * std::powf(10, theexponent - 1);
			thebase = 5;
		}
	}
	else
	{
		if ((log10(cam_w) - theexponent) >= .7)
		{
			spacing = 5 * std::powf(10, theexponent - 1);
			thebase = 5;
		}
		else if ((log10(cam_w) - theexponent) >= .3)
		{
			spacing = 2 * std::powf(10, theexponent - 1);
			thebase = 2;
		}
		else
		{
			spacing = 1 * std::powf(10, theexponent - 1);
			thebase = 1;
		}
	}

	//subgrids
	//smaller gridpoints
	//2s fill in 10, 5s fill in 2, 10s fill in 5
	//make grids where the spacing is from above, then maybe remove the whole one if you want
	float sub_spacing = 0.0f;
	if (thebase == 1)
	{
		sub_spacing = 2 * std::powf(10, theexponent - 2);
	}
	else if (thebase == 5)
	{
		sub_spacing = 1 * std::powf(10, theexponent - 1);
	}
	else if (thebase == 2)
	{
		sub_spacing = 5 * std::powf(10, theexponent - 2);
	}
	//std::cout << "base: " << thebase << "  theexponent: " << theexponent << std::endl;
	//std::cout << sub_spacing << std::endl;
	
	fillgriddata(spacing, grid_points, true, theexponent, sub_spacing);
	writegridbuffer();

	fillgriddata(sub_spacing, subgrid_points, false, 0, 0);
	
	writesubgridbuffer();
}

void Graph2D::fillgriddata(float spacing, std::vector<float>& point_data, bool add_numbers, int exponent, float sub_space)
{
	point_data.clear();
	if (add_numbers)
		rendertextqueue.clear();
	float format_numbers = std::max(exponent - 1, 1);
	int precision_number = 0;
	float extra_scaling = std::max(log10(cam_w) - (int)log10(cam_w), .8f);
	if (exponent <= 0)
	{
		format_numbers = std::abs(exponent - 1);
		precision_number = std::abs(exponent - 1);
		extra_scaling = 1;// std::max(log10(cam_w) - std::floor(log10(cam_w)), .5f);
	}
	float scale = (extra_scaling * spacing) / (150.0f * format_numbers * .6);

	float max_pos = 100000.0f;
	int min_index = ((cam_center.x - cam_w / 2.0f)) / spacing;
	int max_index = ((cam_center.x + cam_w / 2.0f)) / spacing;
	int startingindex = std::min(min_index, max_index);
	int endingindex = std::max(min_index, max_index);
	for (float x = spacing * startingindex; x <= spacing * endingindex; x += spacing)
	{
		point_data.push_back(x);
		point_data.push_back(-max_pos);
		point_data.push_back(0.0f);

		point_data.push_back(x);
		point_data.push_back(max_pos);
		point_data.push_back(0.0f);

		//we need this number to be at (x,0) in the world but also scaled to the same relative size
		if (add_numbers)
		{			
			std::stringstream stream;
			stream << std::fixed << std::setprecision(precision_number) << x;
			std::string s = stream.str();

			//if the xaxis is out of the camera either render them at top or bottom of screen
			if (cam_center.y - cam_h / 2.0f > 0)
			{
				//render at bottom
				rendertextqueue.push_back(GraphicsManager::RenderText(s, (x - (sub_space * 1.0f)) - cam_center.x, (cam_center.y - cam_h / 2.0f) - cam_center.y, scale, glm::vec3(0, 0, 0)));
			}
			else if (cam_center.y + cam_h / 2.0f < 0)
			{
				//render at top
				rendertextqueue.push_back(GraphicsManager::RenderText(s, (x - (sub_space * 1.0f)) - cam_center.x, (cam_center.y - 2*sub_space + cam_h / 2.0f) - cam_center.y, scale, glm::vec3(0, 0, 0)));
			}
			else
			{
				//render on axis
				rendertextqueue.push_back(GraphicsManager::RenderText(s, (x - (sub_space * 1.0f)) - cam_center.x, (0 - (sub_space * 2.0f)) - cam_center.y, scale, glm::vec3(0, 0, 0)));
			}
			
		}
	}

	min_index = ((cam_center.y - cam_h / 2.0f)) / spacing;
	max_index = ((cam_center.y + cam_h / 2.0f)) / spacing;
	startingindex = std::min(min_index, max_index);
	endingindex = std::max(min_index, max_index);
	for (float y = spacing * startingindex; y <= spacing * endingindex; y += spacing)
	{
		point_data.push_back(-max_pos);
		point_data.push_back(y);
		point_data.push_back(0.0f);

		point_data.push_back(max_pos);
		point_data.push_back(y);
		point_data.push_back(0.0f);
		
		if (add_numbers)
		{
			std::stringstream stream;
			stream << std::fixed << std::setprecision(precision_number) << y;
			std::string s = stream.str();

			//if the yaxis is out of the camera either render them at top or bottom of screen
			if (cam_center.x - cam_w / 2.0f > 0)
			{
				//render at left
				rendertextqueue.push_back(GraphicsManager::RenderText(s, (cam_center.x - cam_w/2.0f)-cam_center.x, (y - sub_space / 2.0f) - cam_center.y, scale, glm::vec3(0, 0, 0)));
			}
			else if (cam_center.x + cam_w / 2.0f < 0)
			{
				//render at right
				rendertextqueue.push_back(GraphicsManager::RenderText(s, (cam_center.x + cam_w / 2.0f - (cam_w/15.0f)) - cam_center.x, (y - sub_space / 2.0f) - cam_center.y, scale, glm::vec3(0, 0, 0)));
			}
			else
			{
				//render on axis
				rendertextqueue.push_back(GraphicsManager::RenderText(s, 0-cam_center.x, (y - sub_space/2.0f) - cam_center.y, scale, glm::vec3(0, 0, 0)));
			}
		}
	}
}

void Graph2D::setviewport(float x, float y, float width, float height)
{
	viewport_data = glm::vec4(x, y, width, height);
}


void Graph2D::plotpoint(int line, const std::pair<float, float>& point)
{
	data_points[line].push_back(point.first);
	data_points[line].push_back(point.second);
	data_points[line].push_back(0);

	writegpubuffer(line);
}
/*
void Graph2D::plotpoints(const std::vector<std::pair<float, float>>& points)
{
	for (int i = 0; i < points.size(); i++)
	{
		data_points.push_back(points[i].first);
		data_points.push_back(points[i].second);
		data_points.push_back(0);
	}

	writegpubuffer();
}
*/
void Graph2D::plotpoints(int line, const std::vector<float>& points)
{
	for (int i = 0; i < points.size(); i+=2)
	{
		data_points[line].push_back(points[i]);
		data_points[line].push_back(points[i+1]);
		data_points[line].push_back(0);
	}

	writegpubuffer(line);
}

void Graph2D::clearalllines()
{
	for (int i = 0; i < data_points.size(); i++)
	{
		data_points[i].clear();
		writegpubuffer(i);
	}
}

void Graph2D::clearline(int line)
{
	data_points[line].clear();
}

void Graph2D::createexponential()
{
	auto lambda = [](float x) -> float { return std::exp(x); };
	createfunction(0, 10000, 1, lambda);
}

void Graph2D::createxsquared()
{
	auto lambda = [](float x) -> float { return x*x; };
	createfunction(0, 1000, 1, lambda);
}

void Graph2D::createsin(float A, float w)
{
	auto lambda = [=](float x) -> float { return A*std::sin(w*x); };
	createfunction(0, 10000, 1, lambda);
}

void Graph2D::createlinear()
{
	auto lambda = [](float x) -> float { return x; };
	createfunction(0, 50, 1, lambda);
}

void Graph2D::createfunction(int line, float x_range, float x_precision, std::function<float(float)> graph_function)
{
	data_points[line].clear();
	for (float x = -x_range; x <= x_range; x += x_precision)
	{
		data_points[line].push_back(x);
		data_points[line].push_back(graph_function(x));
		data_points[line].push_back(0);//z is 0
	}

	writegpubuffer(line);
}

void Graph2D::createfunction(int line, float a, float b, float x_precision, std::function<float(float)> graph_function)
{
	data_points[line].clear();
	for (float x = a; x <= b + .0001; x += x_precision)
	{
		data_points[line].push_back(x);
		data_points[line].push_back(graph_function(x));
		data_points[line].push_back(0);//z is 0
	}

	writegpubuffer(line);
}

void Graph2D::DrawToTexture()
{	
	cam_m = glm::translate(glm::mat4(1.0f), glm::vec3(-cam_center.x, -cam_center.y,0));
	proj_m = glm::ortho(-(cam_w / 2.0f), (cam_w / 2.0f), -(cam_h / 2.0f), (cam_h / 2.0f),
		                -1.0f, 1.0f);


	graphics_manager.setFBOstate(fbo_id, background_color, proj_m, cam_m, resolution_w, resolution_h);
	graphics_manager.Begin();

	//subgrids
	GLuint subgrid_count = subgrid_points.size();
	GLuint subgrid_width = 1;
	graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(200,200,200),
										glm::vec4(0.0), {glm::vec2(0,0),glm::vec2(0,0),glm::vec2(0,0)}, 0,
										vao_subgrids,subgrid_width,subgrid_count });
	//grids
	GLuint grid_count = grid_points.size();
	GLuint grid_width = 1.5;
	graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(50,50,50),
										glm::vec4(0.0), {glm::vec2(0,0),glm::vec2(0,0),glm::vec2(0,0)}, 0,
										vao_grids,grid_width,grid_count });

	//axis
	graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(1,0,0), glm::vec4(0,-50000,0,50000) });
	graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(1,0,0), glm::vec4(-50000,0,50000,0) });

	graphics_manager.Render();
	graphics_manager.ClearSubmissions();

	//axis numbers
	for (const auto& text : rendertextqueue)
	{
		graphics_manager.RenderTextFont(text.text, text.x, text.y, text.scale, text.color);
	}
	//graphics_manager.RenderTextFont("50.00", -cam_center.x, -cam_center.y, 1, glm::vec3(0,0,0));
	
	//function
	for (int i = 0; i < active_lines.size(); i++)
	{
		if (!data_points[i].empty() && active_lines[i])
		{
			GLuint vertex_count = data_points[i].size();
			GLuint line_width = 4;
			graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE_STRIP, SHADER_TYPE::MAIN, glm::mat4(1.0f), line_colors[i],
													glm::vec4(0.0), {glm::vec2(0,0),glm::vec2(0,0),glm::vec2(0,0)}, 0,
													vao_points[i],line_width,vertex_count });
		}
	}
	

	//points
	//glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(200, 200, 0));
	//mat = glm::scale(mat, glm::vec3(10, 10, 1));
	//graphics_manager.SubmitDraw({ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
	//mat = glm::translate(glm::mat4(1.0f), glm::vec3(220, 220, 0));
	//mat = glm::scale(mat, glm::vec3(10, 10, 1));
	//graphics_manager.SubmitDraw({ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });


	graphics_manager.Render();
	graphics_manager.ClearSubmissions();

	glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo_id);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo_blit_id);

	glBlitFramebuffer(0, 0, resolution_w, resolution_h, 0, 0, resolution_w, resolution_h, GL_COLOR_BUFFER_BIT, GL_NEAREST);

	graphics_manager.restorFBOstate();
}

void Graph2D::writegpubuffer(int line)
{
	glBindVertexArray(vao_points[line]);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_points[line]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * data_points[line].size(), data_points[line].data(), GL_DYNAMIC_DRAW);
	glBindVertexArray(0);
}

void Graph2D::writegridbuffer()
{
	glBindVertexArray(vao_grids);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_grids);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * grid_points.size(), grid_points.data(), GL_DYNAMIC_DRAW);
	glBindVertexArray(0);
}

void Graph2D::writesubgridbuffer()
{
	glBindVertexArray(vao_subgrids);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_subgrids);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * subgrid_points.size(), subgrid_points.data(), GL_DYNAMIC_DRAW);
	glBindVertexArray(0);
}

uint32_t Graph2D::point_index(uint32_t point)
{
	//3 for every (point-1) points. then you want the next 1 +1, but we start at 0 so -1
	return 3 * (point - 1);
}