#pragma once
#include "../pch.h"
#include "ShaderCompiler.h"
#include "../ThirdParty/stb_image.h"
#include "../ThirdParty/imgui/imgui.h"
#include "../ThirdParty/imgui/imgui_impl_opengl3.h"
#include "../ThirdParty/imgui/imgui_impl_glfw.h"
#include <map>
#include <queue>
#include <ft2build.h>
#include FT_FREETYPE_H

enum class GEOMETRY_TYPE : uint32_t { RECTANGLE, TRIANGLE, CIRCLE, LINE, DYN_TRIANGLE, LINE_STRIP };
enum class SHADER_TYPE : uint32_t { MAIN, TEXTURE };
struct DrawSubmit
{
	GEOMETRY_TYPE geometry;
	SHADER_TYPE shader;
	glm::mat4 mat;
	glm::vec3 color;
	glm::vec4 line_points;
	std::array<glm::vec2,3> tri_points;
	GLuint sampler_id;
	GLuint vao_id; //if you want to pass in your own vao to bind to
	GLuint line_width;
	GLuint vertex_count; //how many vertices we have in vao_id
};
//DrawSubmit{GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(0,0,0), glm::vec4(0.0), {glm::vec2(0,0),glm::vec2(0,0),glm::vec2(0,0)},0,0,0,0}

//handles glfw window as well
class GraphicsManager
{
public:
	GraphicsManager()
	{
		window = nullptr;
		window_height = 0;
		window_width = 0;
		DrawSubmissions.reserve(200);
		PlotData.reserve(5);

		proj_m = glm::mat4(1.0f);
		cam_m = glm::mat4(1.0f);

		clear_color = glm::vec4(0.7, 0.7, 0.7, 1.0);
		default_line_width = 2;
	}

	void CleanUp()
	{
		if(window)
			glfwDestroyWindow(window);

		glfwTerminate();
	}

	int InitWindow(int width, int height, std::string title, std::vector<std::string> icon_image_paths)
	{
		window_width = width;
		window_height = height;

		if (glfwInit() == GLFW_FALSE)
		{
			std::cout << "glfw init failed\n";
			return EXIT_FAILURE;
		}

		glfwSetErrorCallback(error_callback);

		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
		glfwWindowHint(GLFW_SAMPLES, 8);
		glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

		window = glfwCreateWindow(window_width, window_height, title.data(), nullptr, nullptr);
		if (!window)
		{
			std::cout << "window failed\n";
			return EXIT_FAILURE;
		}

		glfwMakeContextCurrent(window);
		glfwSwapInterval(1);

		if (!icon_image_paths.empty())
		{
			GLFWimage images[3];
			int bitDepth = -1;
			images[0].pixels = stbi_load(icon_image_paths[0].data(), &images[0].width, &images[0].height, &bitDepth, 0);
			images[1].pixels = stbi_load(icon_image_paths[1].data(), &images[1].width, &images[1].height, &bitDepth, 0);
			images[2].pixels = stbi_load(icon_image_paths[1].data(), &images[2].width, &images[2].height, &bitDepth, 0);
			if (images[0].pixels == nullptr || images[1].pixels == nullptr || images[2].pixels == nullptr) {
				std::cout << "can't load: " << icon_image_paths[0] << "\n";
			}
			glfwSetWindowIcon(window, 3, images);
			if (images[0].pixels)
				stbi_image_free(images[0].pixels);
			if (images[1].pixels)
				stbi_image_free(images[1].pixels);
			if (images[2].pixels)
				stbi_image_free(images[2].pixels);
		}

		return EXIT_SUCCESS;
	}

	int InitGL()
	{
		if (glewInit() != GLEW_OK)
		{
			std::cout << "glew init failed\n";
			return EXIT_FAILURE;
		}

		//IMGUI
		const char* glsl_version = "#version 460";
		ImGui::CreateContext();
		ImGuiIO& io = ImGui::GetIO(); (void)io;
		ImGui::StyleColorsDark();
		ImGui_ImplGlfw_InitForOpenGL(window, true);
		ImGui_ImplOpenGL3_Init(glsl_version);

		CreateBuffers();
		checkError("createdBuffers");
		CreatePrograms();
		checkError("createdprograms");
		SetPipelineDefault();
		checkError("setpipeline");
		setupFonts();
		checkError("fonts");
		proj_m = glm::ortho(-(window_width / 2.0f), window_width / 2.0f, -(window_height / 2.0f), window_height / 2.0f, -1.0f, 1.0f);

		return EXIT_SUCCESS;
	}

	void Begin()
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}

	void SubmitDraw(DrawSubmit submission)
	{
		DrawSubmissions.push_back(submission);
	}

	void SubmitPlot(std::pair<float, float>* data, uint32_t indices)
	{

	}

	void CreatePlot(float* data, uint32_t indices, glm::mat4 m_plotmatrix)
	{
		std::vector<float> vdata;
		for (int i = 0; i < indices; i++)
		{
			vdata.push_back(data[i]);
			vdata.push_back((float)i);
			vdata.push_back(0.0);
		}

		PlotData = vdata;

		glBindVertexArray(vao_plot);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_plot);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vdata.size(), vdata.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
		glEnableVertexAttribArray(0);

		plotmatrix = m_plotmatrix;
	}
	void DeletePlot(int index)
	{
		//
	}


	void submitPlot(int index)
	{
		plot_submitted = true;
	}

	void Render()
	{
		glm::mat4 vp = proj_m * cam_m;

		for (auto& sub : DrawSubmissions)
		{
			switch (sub.shader)
			{
			case SHADER_TYPE::MAIN: glUseProgram(program_main); break;
			case SHADER_TYPE::TEXTURE: glUseProgram(program_texture); break;
			};
			bool textured = (sub.shader == SHADER_TYPE::TEXTURE);

			if (sub.geometry == GEOMETRY_TYPE::LINE)
			{
				if (sub.vao_id == 0)
				{
					std::vector<float> data = {
						sub.line_points.x, sub.line_points.y, 0.0,
						sub.line_points.z,sub.line_points.w, 0.0
					};

					glBindVertexArray(vao_line);
					glBindBuffer(GL_ARRAY_BUFFER, vbo_line);
					glBufferData(GL_ARRAY_BUFFER, sizeof(float) * data.size(), data.data(), GL_DYNAMIC_DRAW);
				}
				else
				{
					glBindVertexArray(sub.vao_id);
				}

				if (sub.mat != glm::mat4(1.0f)) {
					glUniformMatrix4fv(program_main_mvp, 1, GL_FALSE, glm::value_ptr(vp * sub.mat));
				}
				else
				{
					glUniformMatrix4fv(program_main_mvp, 1, GL_FALSE, glm::value_ptr(vp * glm::mat4(1.0f)));
				}
				glUniform3fv(program_main_col, 1, glm::value_ptr(sub.color));

				if(sub.line_width == 0)
					glLineWidth(default_line_width);
				else
					glLineWidth(sub.line_width);
				
				if(sub.vao_id == 0)
					glDrawArrays(GL_LINES, 0, 2);
				else
					glDrawArrays(GL_LINES, 0, sub.vertex_count / 3);
				
				glLineWidth(default_line_width);
			}
			else if (sub.geometry == GEOMETRY_TYPE::LINE_STRIP)
			{
				glBindVertexArray(sub.vao_id);

				if (sub.mat != glm::mat4(1.0f)) {
					glUniformMatrix4fv(program_main_mvp, 1, GL_FALSE, glm::value_ptr(vp * sub.mat));
				}
				else
				{
					glUniformMatrix4fv(program_main_mvp, 1, GL_FALSE, glm::value_ptr(vp * glm::mat4(1.0f)));
				}
				glUniform3fv(program_main_col, 1, glm::value_ptr(sub.color));
				
				glLineWidth(sub.line_width);
				glDrawArrays(GL_LINE_STRIP, 0, sub.vertex_count / 3);
				glLineWidth(default_line_width);//back to default for now
			}
			else if (sub.geometry == GEOMETRY_TYPE::DYN_TRIANGLE)
			{
				std::vector<float> data = {
					sub.tri_points[0].x, sub.tri_points[0].y, 0.0,
					sub.tri_points[1].x, sub.tri_points[1].y, 0.0,
					sub.tri_points[2].x, sub.tri_points[2].y, 0.0
				};

				glBindVertexArray(vao_dytri);
				glBindBuffer(GL_ARRAY_BUFFER, vbo_dytri);
				glBufferData(GL_ARRAY_BUFFER, sizeof(float) * data.size(), data.data(), GL_DYNAMIC_DRAW);

				glUniformMatrix4fv(program_main_mvp, 1, GL_FALSE, glm::value_ptr(vp*glm::mat4(1.0f)));
				glUniform3fv(program_main_col, 1, glm::value_ptr(sub.color));

				glDrawArrays(GL_TRIANGLES, 0, 3);
			}
			else
			{
				GLsizei index_count = 0;
				switch (sub.geometry)
				{
				case GEOMETRY_TYPE::RECTANGLE: (textured) ? glBindVertexArray(vao_rect_uv) : glBindVertexArray(vao_rect); index_count = 6; break;
				case GEOMETRY_TYPE::TRIANGLE:  glBindVertexArray(vao_tri); index_count = 3;  break;
				case GEOMETRY_TYPE::CIRCLE:  glBindVertexArray(vao_circle); index_count = circle_index_count;  break;
				};

				switch (sub.shader)
				{
				case SHADER_TYPE::MAIN: 
					glUniformMatrix4fv(program_main_mvp, 1, GL_FALSE, glm::value_ptr(vp * sub.mat));
					glUniform3fv(program_main_col, 1, glm::value_ptr(sub.color));
					break;
				case SHADER_TYPE::TEXTURE: 
					glActiveTexture(GL_TEXTURE0);
					glBindTexture(GL_TEXTURE_2D, sub.sampler_id);
					glUniformMatrix4fv(program_texture_mvp, 1, GL_FALSE, glm::value_ptr(vp * sub.mat));
					break;
				};

				glDrawElements(GL_TRIANGLES, index_count, GL_UNSIGNED_INT, 0);
				glBindTexture(GL_TEXTURE_2D, 0);
			}
		}

		if (plot_submitted && !PlotData.empty())
		{
			glUseProgram(program_main);

			glBindVertexArray(vao_plot);

			glUniformMatrix4fv(program_main_mvp, 1, GL_FALSE, glm::value_ptr(vp * plotmatrix));
			glUniform3fv(program_main_col, 1, glm::value_ptr(glm::vec3(0, 0, 0)));
			
			glPointSize(5);
			glDrawArrays(GL_POINTS, 0, PlotData.size()-1);
		}

		checkError("finished render");
	}
	
	void RenderTextFont(std::string text, float x, float y, float scale, glm::vec3 color)
	{
		glUseProgram(program_font);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glActiveTexture(GL_TEXTURE0);
		glBindVertexArray(VAO_font);

		glm::mat4 p_mat = proj_m;
		glUniformMatrix4fv(program_font_mvp, 1, GL_FALSE, glm::value_ptr(p_mat));

		glUniform3f(program_font_col, color.x, color.y, color.z);

		// iterate through all characters
		std::string::const_iterator c;
		for (c = text.begin(); c != text.end(); c++)
		{
			Character ch = Characters[*c];

			float xpos = x + ch.Bearing.x * scale;
			float ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

			float w = ch.Size.x * scale;
			float h = ch.Size.y * scale;
			// update VBO for each character
			float vertices[6][4] = {
				{ xpos,     ypos + h,   0.0f, 0.0f },
				{ xpos,     ypos,       0.0f, 1.0f },
				{ xpos + w, ypos,       1.0f, 1.0f },

				{ xpos,     ypos + h,   0.0f, 0.0f },
				{ xpos + w, ypos,       1.0f, 1.0f },
				{ xpos + w, ypos + h,   1.0f, 0.0f }
			};
			// render glyph texture over quad
			glBindTexture(GL_TEXTURE_2D, ch.TextureID);
			// update content of VBO memory
			glBindBuffer(GL_ARRAY_BUFFER, VBO_font);
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			// render quad
			glDrawArrays(GL_TRIANGLES, 0, 6);
			// now advance cursors for next glyph (note that advance is number of 1/64 pixels)
			x += (ch.Advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64)
		}
		glBindVertexArray(0);
		glBindTexture(GL_TEXTURE_2D, 0);
		glDisable(GL_BLEND);
	}

	void End()
	{
		glfwSwapBuffers(window);
		DrawSubmissions.clear();
		plot_submitted = false;
	}

	void ClearSubmissions()
	{
		DrawSubmissions.clear();
	}

	void setFBOstate(GLuint fbo_id, glm::vec4 clear_col, glm::mat4 proj, glm::mat4 cam, int w, int h)
	{
		glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);
		glClearColor(clear_col.x, clear_col.y, clear_col.z, clear_col.w);

		proj_m_save = proj_m;
		proj_m = proj;
		cam_m_save = cam_m;
		cam_m = cam;

		glViewport(0, 0, w, h);
		default_line_width = 3;
		glLineWidth(default_line_width);
	}

	void restorFBOstate()
	{
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);

		proj_m = proj_m_save;
		cam_m = cam_m_save;
		glViewport(0, 0, window_width, window_height);
		default_line_width = 2;
		glLineWidth(2);
	}

	void WireFrame(bool enable)
	{
		(enable) ? glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) : glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		window_width = width;
		window_height = height;
		proj_m = glm::ortho(-(window_width / 2.0f), window_width / 2.0f, -(window_height / 2.0f), window_height / 2.0f, -1.0f, 1.0f);
		proj_m_save = proj_m;

		glViewport(0, 0, width, height);
	}

	GLFWwindow* getWindow() { return window; }

	uint32_t getWidth() { return window_width; }
	uint32_t getHeight() { return window_height; }

	//converts the window space to clip (-1,1) space
	glm::vec2 WindowToClip(glm::vec2 pos, bool cursor_used)
	{
		glm::vec2 clip = pos;
		clip.x /= window_width;//(0,1)
		clip.x = clip.x * 2.0 - 1.0; //(-1,1)

		clip.y /= window_height;
		clip.y = clip.y * 2.0 - 1.0;
		if(cursor_used)
			clip.y *= -1.0;

		return clip;
	}

	//converts the (0,0) at top left (W,H) bottom right, to the origin being at the center
	glm::vec2 WindowToWorld(glm::vec2 pos)
	{
		//what you see in your window is a consequence of the projection matrix pretty much (technically all the transforsm in the whole pipeline
		//convert your window to ndc, then scale that by the ortho_box size

		glm::vec2 ndc = WindowToClip(pos, true);

		glm::vec2 worldc;
		worldc.x = ndc.x * (window_width / 2.0f);
		worldc.y = ndc.y * (window_height / 2.0f);

		return worldc;
	}

	glm::vec2 WindowToClipGeneral(glm::vec2 pos, int width, int height, bool cursor_used)
	{
		//technically the center is the viewport center, and the width/height are the viewport width height
		//you need to find where your mouse is not in terms of the window but the viewport

		glm::vec2 clip = pos;
		clip.x /= width;//(0,1)
		clip.x = clip.x * 2.0 - 1.0; //(-1,1)

		clip.y /= height;
		clip.y = clip.y * 2.0 - 1.0;
		if (cursor_used)
			clip.y *= -1.0;

		return clip;
	}

	glm::vec2 WindowToWorldGeneral(glm::vec2 pos, int window_w, int window_h, glm::vec4 viewport, float cam_w, float cam_h, float cam_x, float cam_y)
	{
		//we need to convert to NDC
		glm::vec2 ndc = WindowToClipGeneral(pos, window_width, window_h, true);

		//now we need to translate and scale to the viewport (which lives in ndc)
		glm::vec2 viewport_pos = glm::vec2(viewport.x, viewport.y) / glm::vec2(window_w/2.0f, window_h/2.0f);
		glm::vec2 viewport_dim = glm::vec2(viewport.z,viewport.w) / glm::vec2(window_w, window_h);
		
		//move it to new origin, if its to the right of the viewport origina it has positive x position
		ndc = ndc - viewport_pos;
	
		//scale it, if its half the size of window it gets scaled double
		ndc /= viewport_dim;

		//then we need to convert to camera space, the camera has a width and height 
		glm::vec2 worldc;
		worldc.x = ndc.x * (cam_w / 2.0f);
		worldc.y = ndc.y * (cam_h / 2.0f);

		//the camera is someone in the world so you need to offset that
		worldc.x += cam_x;
		worldc.y += cam_y;

		return worldc;
	}

	void setProjM(glm::mat4 m)
	{
		proj_m = m;
	}

	struct RenderText {
		std::string text;
		float x;
		float y;
		float scale;
		glm::vec3 color;

		RenderText(std::string _text, float _x, float _y, float _scale, glm::vec3 _color) :
			text(_text), x(_x), y(_y), scale(_scale), color(_color) {}
	};
private:
	GLFWwindow* window;
	uint32_t window_width, window_height;

	glm::mat4 proj_m;
	glm::mat4 proj_m_save;
	glm::mat4 cam_m;
	glm::mat4 cam_m_save;
	//draw info:shape, shader, inputs, state
	//to optimize I can sort these by shaders and geometry etc, right now its fine
	std::vector<DrawSubmit> DrawSubmissions;
	std::vector<float> PlotData;
	bool plot_submitted;
	glm::mat4 plotmatrix;

	glm::vec4 clear_color;
	GLuint default_line_width;

	GLuint vao_rect, vbo_rect, ibo_rect;
	GLuint vao_rect_uv, vbo_rect_uv, ibo_rect_uv;
	GLuint vao_tri, vbo_tri, ibo_tri;
	GLuint vao_circle, vbo_circle, ibo_circle;
	GLuint vao_line, vbo_line;
	GLuint vao_dytri, vbo_dytri;
	GLuint circle_index_count;
	GLuint vao_plot, vbo_plot;

	//fonts
	struct Character {
		unsigned int TextureID;  // ID handle of the glyph texture
		glm::ivec2   Size;       // Size of glyph
		glm::ivec2   Bearing;    // Offset from baseline to left/top of glyph
		unsigned int Advance;    // Offset to advance to next glyph
	};

	unsigned int VAO_font, VBO_font;
	std::map<char, Character> Characters;
	FT_Library ft;
	FT_Face face;
	std::queue<RenderText> rendertextqueue;
	GLuint program_font, program_font_mvp, program_font_col;

	GLuint program_main, program_main_mvp, program_main_col;
	GLuint program_texture, program_texture_mvp;
private:
	
	void CreateBuffers()
	{
		//RECT
		std::vector<float> rect_pos =
		{
			-.5,.5,0,
			-.5,-.5,0,
			.5,-.5,0,
			.5,.5,0
		};

		std::vector<uint32_t> rect_idx =
		{
			3,0,1,
			3,1,2
		};

		glGenVertexArrays(1, &vao_rect);
		glBindVertexArray(vao_rect);
		glGenBuffers(1, &vbo_rect);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_rect);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * rect_pos.size(), rect_pos.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
		glEnableVertexAttribArray(0);
		//glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(sizeof(float) * 3)); //uv attribute
		//glEnableVertexAttribArray(1);

		glGenBuffers(1, &ibo_rect);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_rect);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * rect_idx.size(), rect_idx.data(), GL_STATIC_DRAW);

		//RECT_UV
		std::vector<float> rect_uv =
		{
			-.5,.5,0,  0, 1,
			-.5,-.5,0, 0, 0,
			.5,-.5,0,  1, 0,
			.5,.5,0,   1, 1
		};

		glGenVertexArrays(1, &vao_rect_uv);
		glBindVertexArray(vao_rect_uv);
		glGenBuffers(1, &vbo_rect_uv);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_rect_uv);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * rect_uv.size(), rect_uv.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), 0); //position attribute
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(sizeof(float) * 3)); //uv attribute
		glEnableVertexAttribArray(1);

		glGenBuffers(1, &ibo_rect_uv);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_rect_uv);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * rect_idx.size(), rect_idx.data(), GL_STATIC_DRAW);

		//TRI

		float a = 1;
		float height = a * sin(glm::radians(60.0f));	
		float width = a/2.0f;
		std::vector<float> tri_pos =
		{
			0,height/2.0f,0,
			-width,-height/3.0f,0,
			width,-height/3.0f,0
		};

		std::vector<uint32_t> tri_idx =
		{
			0,1,2
		};

		glGenVertexArrays(1, &vao_tri);
		glBindVertexArray(vao_tri);
		glGenBuffers(1, &vbo_tri);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_tri);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * tri_pos.size(), tri_pos.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
		glEnableVertexAttribArray(0);
		//glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(sizeof(float) * 3)); //uv attribute
		//glEnableVertexAttribArray(1);

		glGenBuffers(1, &ibo_tri);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_tri);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * tri_idx.size(), tri_idx.data(), GL_STATIC_DRAW);

		//CIRCLE
		std::vector<float> circle_pos;
		std::vector<uint32_t> circle_idx;

		float N = 30;
		float increment = 360 / N;

		circle_pos.push_back(0);
		circle_pos.push_back(0);
		circle_pos.push_back(0);
		for (float theta = 0; theta < 360; theta += increment)
		{
			circle_pos.push_back(cos(glm::radians(theta)));
			circle_pos.push_back(sin(glm::radians(theta)));
			circle_pos.push_back(0);
		}

		for (int i = 1; i < N; i++)
		{
			circle_idx.push_back(i);
			circle_idx.push_back(i + 1);
			circle_idx.push_back(0);
		}
		circle_idx.push_back(N);
		circle_idx.push_back(1);
		circle_idx.push_back(0);

		circle_index_count = circle_idx.size();

		glGenVertexArrays(1, &vao_circle);
		glBindVertexArray(vao_circle);
		glGenBuffers(1, &vbo_circle);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_circle);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * circle_pos.size(), circle_pos.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
		glEnableVertexAttribArray(0);
		//glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(sizeof(float) * 3)); //uv attribute
		//glEnableVertexAttribArray(1);

		glGenBuffers(1, &ibo_circle);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_circle);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * circle_idx.size(), circle_idx.data(), GL_STATIC_DRAW);

		//LINE
		glGenVertexArrays(1, &vao_line);
		glBindVertexArray(vao_line);
		glGenBuffers(1, &vbo_line);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_line);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
		glEnableVertexAttribArray(0);

		//DYTRI
		glGenVertexArrays(1, &vao_dytri);
		glBindVertexArray(vao_dytri);
		glGenBuffers(1, &vbo_dytri);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_dytri);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
		glEnableVertexAttribArray(0);

		//PLOT
		//DYTRI
		glGenVertexArrays(1, &vao_plot);
		glBindVertexArray(vao_plot);
		glGenBuffers(1, &vbo_plot);
		glBindBuffer(GL_ARRAY_BUFFER, vbo_plot);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0); //position attribute
		glEnableVertexAttribArray(0);
	}

	void CreatePrograms()
	{
		//main program
		ShaderCompiler::ReadFile("Shaders/main.vert");
		ShaderCompiler::CompileGLShader("Shaders/main.vert", "Shaders/main.frag", program_main);
		glUseProgram(program_main);
		program_main_mvp = glGetUniformLocation(program_main, "mvp");
		program_main_col = glGetUniformLocation(program_main, "col");
		glUseProgram(0);

		//texture program
		ShaderCompiler::ReadFile("Shaders/texture.vert");
		ShaderCompiler::CompileGLShader("Shaders/texture.vert", "Shaders/texture.frag", program_texture);
		glUseProgram(program_texture);
		program_texture_mvp = glGetUniformLocation(program_texture, "mvp");
		glUseProgram(0);

		//font program
		ShaderCompiler::ReadFile("Shaders/font.vert");
		ShaderCompiler::CompileGLShader("Shaders/font.vert", "Shaders/font.frag", program_font);
		glUseProgram(program_font);
		program_font_mvp = glGetUniformLocation(program_font, "projection");
		program_font_col = glGetUniformLocation(program_font, "textColor");
		glUseProgram(0);
	}

	void SetPipelineDefault()
	{
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_SCISSOR_TEST);
		glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
		glEnable(GL_PROGRAM_POINT_SIZE);
		//glEnable(GL_LINE_SMOOTH);

		glDepthFunc(GL_LESS);
		glDepthMask(GL_FALSE);
		glClearDepth(1.0);

		glEnable(GL_MULTISAMPLE);
		//glEnable(GL_BLEND);
		//glBlendEquation(GL_FUNC_ADD);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}

	void setupFonts()
	{
		//FREETYPE
		if (FT_Init_FreeType(&ft))
		{
			std::cout << "Could not init FreeType error\n";
			return;
		}
		FT_Int maj, min, patch;
		FT_Library_Version(ft, &maj, &min, &patch);
		std::cout << "FreeType Version: " << maj << " " << min << " " << patch << std::endl;

		if (FT_New_Face(ft, std::string("Game_Data/Fonts/OpenSans-Regular.ttf").c_str(), 0, &face))
		{
			std::cout << "failed to load font OpenSans-Regular.ttf\n";
			return;
		}
		FT_Set_Pixel_Sizes(face, 0, 48);
		if (FT_Load_Char(face, 'X', FT_LOAD_RENDER))
		{
			std::cout << "failed to load glyph\n";
			return;
		}

		glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // disable byte-alignment restriction

		for (unsigned char c = 0; c < 128; c++)
		{
			// load character glyph 
			if (FT_Load_Char(face, c, FT_LOAD_RENDER))
			{
				std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
				continue;
			}
			// generate texture
			unsigned int texture;
			glGenTextures(1, &texture);
			glBindTexture(GL_TEXTURE_2D, texture);
			glTexImage2D(
				GL_TEXTURE_2D,
				0,
				GL_RED,
				face->glyph->bitmap.width,
				face->glyph->bitmap.rows,
				0,
				GL_RED,
				GL_UNSIGNED_BYTE,
				face->glyph->bitmap.buffer
			);

			// set texture options 
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE_EXT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE_EXT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			// now store character for later use
			Character character = {
				texture,
				glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
				glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
				face->glyph->advance.x
			};
			Characters.insert(std::pair<char, Character>(c, character));
		}
		checkError("created font textures");

		glGenVertexArrays(1, &VAO_font);
		glGenBuffers(1, &VBO_font);
		glBindVertexArray(VAO_font);
		glBindBuffer(GL_ARRAY_BUFFER, VBO_font);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);

		checkError("created font vao");
	}
	public:
	bool checkError(std::string msg)
	{
		GLenum errorCode;
		bool any_error = false;
		while ((errorCode = glGetError()) != GL_NO_ERROR)
		{
			any_error = true;
			std::string error = "";
			switch (errorCode)
			{
			case GL_INVALID_ENUM:                  error = "INVALID_ENUM"; break;
			case GL_INVALID_VALUE:                 error = "INVALID_VALUE"; break;
			case GL_INVALID_OPERATION:             error = "INVALID_OPERATION"; break;
			case GL_STACK_OVERFLOW:                error = "STACK_OVERFLOW"; break;
			case GL_STACK_UNDERFLOW:               error = "STACK_UNDERFLOW"; break;
			case GL_OUT_OF_MEMORY:                 error = "OUT_OF_MEMORY"; break;
			case GL_INVALID_FRAMEBUFFER_OPERATION: error = "INVALID_FRAMEBUFFER_OPERATION"; break;
			}
			std::cout << "GLERROR (" << msg << "): " << error << " " << errorCode << std::endl;
		}
		return any_error;
	}

	static void error_callback(int error_code, const char* description)
	{
		std::cout << "GLFW ERROR " << error_code << ": " << description << std::endl;
	}
};


/* IMGUI TEMPLATE

ImGui_ImplOpenGL3_NewFrame();
ImGui_ImplGlfw_NewFrame();

ImGui::NewFrame();
ImGui::Begin("Polygons");
//ImGui::Checkbox("Show Skeleton", &showSkeleton);
//ImGui::Checkbox("Show Mesh", &showMesh);
//ImGui::Checkbox("Show Animations", &showAnimations);
//ImGui::SliderFloat("Animation Time", &animationTime, 0.0f, 100.0f);
//const char* items2[] = { "LINEAR", "SLERP", "SQUAD" };
//ImGui::Combo("Interpolation Mode", &interpolationMethod, items2, 3);
//const char* items3[] = { "Vampire" };
//ImGui::Combo("Skeletal Model", &current_model, items3, 1);
//const char* items4[] = { "Dancing" };
//ImGui::Combo("Animation Track", &current_animation, items4, 1);
ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
ImGui::End();

ImGui::Render();
ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

*/