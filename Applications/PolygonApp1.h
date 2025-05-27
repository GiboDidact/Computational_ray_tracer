#pragma once
#include "../Graphics/GraphicsManager.h"
#include "../Geometry/Polygon.h"
#include "../Geometry/ArtGalleryTheorem.h"

class PolygonApp1
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		PolygonApp1* handler = reinterpret_cast<PolygonApp1*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<PolygonApp1*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<PolygonApp1*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<PolygonApp1*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<PolygonApp1*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons = {
			"Game_Data/AlO16.png",
			"Game_Data/AlO32.png",
			"Game_Data/AlO48.png"
		};
		if (graphics_manager.InitWindow(1500, 1000, "Polygon Engine", icons) == EXIT_FAILURE)
			return;
		if (graphics_manager.InitGL() == EXIT_FAILURE)
			return;

		glfwSetWindowUserPointer(graphics_manager.getWindow(), this);
		glfwSetFramebufferSizeCallback(graphics_manager.getWindow(), framebuffer_size_callback);
		glfwSetKeyCallback(graphics_manager.getWindow(), key_callback);
		glfwSetCursorPosCallback(graphics_manager.getWindow(), cursor_position_callback);
		glfwSetMouseButtonCallback(graphics_manager.getWindow(), mouse_button_callback);
		glfwSetScrollCallback(graphics_manager.getWindow(), scroll_callback);

		MainLoop();
		ShutDown();
	}

	void MainLoop()
	{
		std::vector<Polygon> polygons;

		bool create_polygon = false;
		int vertex_count = 3;
		int selected_poly = -1;
		bool clicked_off = false;
		int selected_vertex = -1;
		
		glm::vec2 last_mouse_pos(0, 0);
		glm::vec2 middle_mouse_pos(0, 0);
		bool moving = false;
		bool moving_vertex = false;
		bool moving_scale = false;
		bool highlight_ear = false;

		bool draw_triangulation = true;
		bool draw_ears = false;
		bool draw_diagnol = false;
		bool draw_mouth = false;
		bool draw_vertex = false;
		bool view_files = false;

		std::vector<std::string> file_names = Helper::LoadPolygonDirectoryFileNames();
		int load_polygon_id = 0;

		Polygon preview_polygon;
		int generate_polygon_N = 6;

		bool imgui_hovered = false;

		//art gallery
		bool place_guards = false;
		bool guards_begin = true;
		std::vector<glm::vec2> guards;
		std::vector<Polygon> visibility_polys;

		//line intersection test
		bool line_intersection_test = false;
		glm::vec2 a1(-90, 300);
		glm::vec2 b1(-90, 100);
		glm::vec2 a2(-100, -100);
		glm::vec2 b2(0, 0);
		int button_id = 0;
		while (!glfwWindowShouldClose(graphics_manager.getWindow()))
		{
			//INPUT
			input_manager.HandleInput();
			glm::vec2 mouse_pos = graphics_manager.WindowToWorld(input_manager.getMousePosition());
			
			if (line_intersection_test)
			{
				if (input_manager.getMouseButtonOnce().x) {
					button_id = (button_id + 1) % 4;
				}
				if(button_id == 0)
					a1 = mouse_pos;
				else if (button_id == 1)
					b1 = mouse_pos;
				else if (button_id == 2)
					a2 = mouse_pos;
				else if (button_id == 3)
					b2 = mouse_pos;
			}

			//rotate selected poly
			if (selected_poly != -1 && !place_guards)
			{
				if (input_manager.MouseScrollDirection() == 1)
				{
					polygons[selected_poly].Rotate(10);
				}
				else if (input_manager.MouseScrollDirection() == -1)
				{
					polygons[selected_poly].Rotate(-10);
				}
			}

			//scale selected poly
			if (moving_scale && input_manager.getMouseMiddleCurrent() && !place_guards)
			{
				//translate the poly based off offset of original click
				glm::vec2 center = polygons[selected_poly].getCenter();
				float r1 = glm::distance(center, mouse_pos);
				float r2 = glm::distance(center, last_mouse_pos);
				last_mouse_pos = mouse_pos;
				float amount = .01;
				if (r1 < r2)
					polygons[selected_poly].Scale(1.0f - amount);
				else if(r1 > r2)
					polygons[selected_poly].Scale(1.0f + amount);
			}
			else
			{
				moving_scale = false;
			}

			//translate selected poly
			if (moving && input_manager.getMouseButtonCurrent().x)
			{
				glm::vec2 offset = mouse_pos - last_mouse_pos;
				last_mouse_pos = mouse_pos;
				polygons[selected_poly].Shift(offset);
			}
			else
			{
				moving = false;
			}

			//translate selected vertex
			if (moving_vertex && input_manager.getMouseButtonCurrent().x)
			{
				glm::vec2 offset = mouse_pos - last_mouse_pos;
				last_mouse_pos = mouse_pos;
				polygons[selected_poly].ShiftVertex(offset, selected_vertex);
				polygons[selected_poly].TriangulateDiagnolSplitting();
			}
			else
			{ 
				moving_vertex = false;
			}

			if (place_guards)
			{
				if (guards_begin)
				{
					guards_begin = false;
				
					//create guard
					guards.push_back(mouse_pos);

					//create visibility polygon
					visibility_polys.push_back(Polygon());
				}

				guards[guards.size() - 1] = mouse_pos;
				//generate visibility polygon for this guard
				visibility_polys[guards.size() - 1] = ArtGalleryTheorem::GenerateVisibilityPolygon(polygons[selected_poly], guards[guards.size() - 1]);
			}
			else
			{
				if (!guards_begin)
				{
					guards.clear();
					visibility_polys.clear();
				}
				guards_begin = true;
			}

			//left click
			if (imgui_hovered == false && input_manager.getMouseButtonOnce().x)
			{
				if (place_guards)
				{
					//place guard at position 
					guards.push_back(mouse_pos);

					visibility_polys.push_back(Polygon());
				}
				else
				{
					bool vertex_found = false;
					//see if we hit a vertex to move
					if (selected_poly != -1)
					{
						selected_vertex = -1;
						int id = polygons[selected_poly].hitVertex(mouse_pos, 10);
						if (id != -1)
						{
							selected_vertex = id;
							vertex_found = true;
							last_mouse_pos = mouse_pos;
							moving_vertex = true;
						}

						//see if we hit edge instead
						if (id == -1)
						{
							int edgeid = polygons[selected_poly].hitBoundaryNoVertex(mouse_pos, 4);
							if (edgeid != -1)
							{
								polygons[selected_poly].splitEdge(edgeid);
								polygons[selected_poly].TriangulateDiagnolSplitting();
							}
						}
					}

					//see if we find a polygon
					if (vertex_found == false)
					{
						selected_poly = -1;
						for (int i = 0; i < polygons.size(); i++)
						{
							if (polygons[i].IsInside(mouse_pos))
							{
								selected_poly = i;
								break;
							}
						}

						if (selected_poly != -1)
						{
							clicked_off = true;
							last_mouse_pos = mouse_pos;
							moving = true;
						}
						else
						{
							clicked_off = false;
						}
					}
				}
			}

			//right click
			if (input_manager.getMouseButtonOnce().y)
			{
				if (place_guards)
				{
					//remove guard at that position
					for (int i = 0; i < guards.size()-1; i++)
					{
						if (abs(glm::distance(mouse_pos, guards[i])) <= 10)
						{
							guards.erase(guards.begin() + i);
							visibility_polys.erase(visibility_polys.begin() + i);
							break;
						}
					}
				}
				else
				{
					//if we are selecting polygon and right click either destroy vertex or the whole polygon
					if (selected_poly != -1)
					{
						int id = polygons[selected_poly].hitVertex(mouse_pos, 20);
						if (id != -1)
						{
							polygons[selected_poly].joinEdge(id);
							polygons[selected_poly].TriangulateDiagnolSplitting();
						}
						else
						{
							polygons.erase(polygons.begin() + selected_poly);
							selected_poly = -1;
							selected_vertex = -1;
							moving = false;
							moving_vertex = false;
							moving_scale = false;
						}
					}
				}
			}

			//scaling begin
			if (input_manager.getMouseMiddleDown())
			{
				if (selected_poly != -1)
				{
					moving_scale = true;
					last_mouse_pos = mouse_pos;
					middle_mouse_pos = mouse_pos;
				}
			}

			//polygon creation if we hit button
			if (create_polygon)
			{
				create_polygon = false;

				//create a convex polygon with vertex_count vertices
				std::vector<glm::vec2> data;
				float r = 5*vertex_count;
				for (float i = 0; i < 360; i += (360 / vertex_count))
				{
					data.push_back(glm::vec2(r*cos(glm::radians(i)), r*sin(glm::radians(i))));
				}
				polygons.push_back(Polygon(data));
				polygons[polygons.size() - 1].TriangulateDiagnolSplitting(true);
			}

			//DRAW
			graphics_manager.Begin();

			if (line_intersection_test)
			{
				glm::vec3 results = Segment::LineToLineIntersectionPoint2(a1, b1, a2, b2);
				if (results.x != -1)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(results.y, results.z, 0));
					mat = glm::scale(mat, glm::vec3(5, 5, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0)});

					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(255,0,0), glm::vec4(a1.x,a1.y,b1.x,b1.y) });
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(0,0,0), glm::vec4(a2.x,a2.y,b2.x,b2.y) });
				}
				else
				{
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(0,0,0), glm::vec4(a1.x,a1.y,b1.x,b1.y) });
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(0,0,0), glm::vec4(a2.x,a2.y,b2.x,b2.y) });
				}

			}

			for (int i = 0; i < polygons.size(); i++)
			{
				if(i == selected_poly)
					polygons[i].DrawPolygon(graphics_manager, draw_triangulation, draw_vertex, glm::vec3(200,100,0));
				else
					polygons[i].DrawPolygon(graphics_manager, draw_triangulation, draw_vertex);

				if(draw_ears)
					polygons[i].HighlightEars(graphics_manager);
				if(draw_diagnol)
					polygons[i].HighlightDiagnol(graphics_manager);
				if (draw_mouth)
					polygons[i].HighlightMouth(graphics_manager);
			}
			if (place_guards || true)
			{
				//draw visibility polygons
				for (int i = 0; i < visibility_polys.size(); i++)
				{
					visibility_polys[i].DrawPolygon(graphics_manager, false, false, glm::vec3(200, 200, 0), false, false);
				}

				//draw guards
				for (int i = 0; i < guards.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(guards[i].x, guards[i].y, 0));
					mat = glm::scale(mat, glm::vec3(15, 15, 1));
					glm::vec3 col = glm::vec3(50, 150, 150);
					if (i == guards.size() - 1)
						col = glm::vec3(10, 50, 50);
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::TRIANGLE, SHADER_TYPE::MAIN, mat, col });
				}
			}

			if (view_files)
			{
				if (file_names.size() != 0)
				{
					std::string file_name = "poly" + std::to_string(load_polygon_id) + ".txt";
					preview_polygon = Polygon::ImportPolygon(file_name);
					preview_polygon.TriangulateDiagnolSplitting(true);
					preview_polygon.DrawPolygon(graphics_manager, false, false, glm::vec3(50, 50, 50));
				}
			}
			graphics_manager.Render();
			
			//IMGUI
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();
			ImGui::NewFrame();
			//bool lol;
			//ImGui::ShowDemoWindow(&lol);
			
			ImGui::Begin("Polygons");
			if (ImGui::Button("Create Polygon"))
				create_polygon = true;
			ImGui::SliderInt("Vertex Count", &vertex_count, 3, 300);

			ImGui::Checkbox("Draw_Triangulations", &draw_triangulation);
			ImGui::Checkbox("Draw_Ears", &draw_ears);
			ImGui::Checkbox("Draw_Diagnols", &draw_diagnol);
			ImGui::Checkbox("Draw_Mouths", &draw_mouth);
			ImGui::Checkbox("Show First Two Vertex (Red-->Green-->Blue)", &draw_vertex);

			if (ImGui::Button("Save Selected Polygon"))
			{
				if (polygons.size() > 0 && selected_poly != -1)
				{
					std::string file_name = "poly" + std::to_string(file_names.size()) + ".txt";
					polygons[selected_poly].ExportPolygon(file_name);

					file_names = Helper::LoadPolygonDirectoryFileNames();
				}
			}
			imgui_hovered = ImGui::IsWindowHovered();

			if (ImGui::Button("Load Polygon"))
			{
				if (file_names.size() != 0)
				{
					std::string file_name = "poly" + std::to_string(load_polygon_id) + ".txt";
					Polygon p = Polygon::ImportPolygon(file_name);

					p.TriangulateDiagnolSplitting(true);
					polygons.push_back(p);

					view_files = false;
				}
			}
			ImGui::SliderInt("File to load", &load_polygon_id, 0, file_names.size() - 1);

			ImGui::Checkbox("View Polygons Files", &view_files);

			if (ImGui::Button("Generate Random Polygon"))
			{
				Polygon p = Polygon::GenerateRandomPolygon(generate_polygon_N);
				p.TriangulateDiagnolSplitting(true);
				polygons.push_back(p);
			}
			ImGui::SliderInt("Generate Vertex Count", &generate_polygon_N, 3, 300);

			if (selected_poly != -1)
			{
				ImGui::Checkbox("Place Guards", &place_guards);

				if (ImGui::Button("Generate Guards(3-Color)"))
				{
					guards.clear();
					visibility_polys.clear();
					guards = ArtGalleryTheorem::PlaceGuardsColorTheorem(polygons[selected_poly]);
					std::cout << guards.size() << std::endl;
					for (int i = 0; i < guards.size(); i++)
					{
						visibility_polys.push_back(ArtGalleryTheorem::GenerateVisibilityPolygon(polygons[selected_poly], guards[i]));
					}
				}
				if (ImGui::Button("Generate Guards(Every Triangle"))
				{
					guards.clear();
					visibility_polys.clear();
					guards = ArtGalleryTheorem::PlaceGuardsTriangles(polygons[selected_poly]);
					for (int i = 0; i < guards.size(); i++)
					{
						visibility_polys.push_back(ArtGalleryTheorem::GenerateVisibilityPolygon(polygons[selected_poly], guards[i]));
					}
					std::cout << guards.size() << std::endl;
				}
				if (ImGui::Button("Generate Guards(Every Other Triangle"))
				{
					guards.clear();
					visibility_polys.clear();
					guards = ArtGalleryTheorem::PlaceGuardsTrianglesShared(polygons[selected_poly]);
					for (int i = 0; i < guards.size(); i++)
					{
						visibility_polys.push_back(ArtGalleryTheorem::GenerateVisibilityPolygon(polygons[selected_poly], guards[i]));
					}
					std::cout << guards.size() << std::endl;
				}
			}
			else
			{
				guards.clear();
				visibility_polys.clear();
			}

			//polygon info
			for (int i = 0; i < polygons.size(); i++)
			{
				std::string text1 = "polygon " + std::to_string(i);
				ImGuiTreeNodeFlags flag = 0;
				if (i == selected_poly) {
					ImGui::SetNextTreeNodeOpen(true);
					flag = ImGuiTreeNodeFlags_Selected;
				}
				else if (clicked_off)
				{
					ImGui::SetNextTreeNodeOpen(false);
				}

				if (ImGui::TreeNodeEx(text1.data(), flag))
				{
					std::string text1 = "vertex count: " + std::to_string(polygons[i].getvertexcount());
					text1 += " triangle count: " + std::to_string(polygons[i].gettrianglecount());
					ImGui::Text(text1.data());

					std::string text2 = "area: " + std::to_string(polygons[i].CalculateArea());
					ImGui::Text(text2.data());

					text2 = "internal angles: " + std::to_string(polygons[i].CalculateInternalAngles());
					ImGui::Text(text2.data());

					text2 = "isvalid: " + std::to_string(polygons[i].isValidPolygon());
					ImGui::Text(text2.data());

					text2 = "IsConvex: " + std::to_string(polygons[i].IsConvex());
					ImGui::Text(text2.data());

					text2 = "reflex vertices: " + std::to_string(polygons[i].getreflexcount());
					ImGui::Text(text2.data());

					text2 = "IsCCW: " + std::to_string(polygons[i].CalculateCCW());
					ImGui::Text(text2.data());

					if (ImGui::Button("Triangulate"))
						polygons[i].TriangulateDiagnolSplitting(true, true);
					if (ImGui::Button("Calculate Ears"))
						polygons[i].CalculateEars(true);
					if (ImGui::Button("Calculate Diagnol"))
						polygons[i].FindDiagnol();
					if (ImGui::Button("Calculate Mouth"))
						polygons[i].CalculateMouth();
					if (ImGui::Button("Make CCW"))
						polygons[i].MakeCCW();

					ImGui::TreePop();
				}
			}

			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
			ImGui::End();
			
			ImGui::Render();
			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

			graphics_manager.End();
		}
	}

	void ShutDown()
	{
		graphics_manager.CleanUp();
	}

	InputManager& GetInputManager() { return input_manager; }
	GraphicsManager& GetGraphicsManager() { return graphics_manager; }

private:
	GraphicsManager graphics_manager;
	InputManager input_manager;
};