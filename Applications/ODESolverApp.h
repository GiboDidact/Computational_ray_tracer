#pragma once
#include "../Graphics/GraphicsManager.h"
#include "../Graphics/Graph2D.h"

class ODESolverApp
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		ODESolverApp* handler = reinterpret_cast<ODESolverApp*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<ODESolverApp*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<ODESolverApp*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<ODESolverApp*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<ODESolverApp*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons;
		if (graphics_manager.InitWindow(1850, 1850 / 1.5f, "ODESolver App", icons) == EXIT_FAILURE)
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
	 
	std::vector<float> eulersmethod(float delta_x, float x_0, float y_0, float range, std::function<float(float,float)> y_derivative)
	{
		//delta_x is precision for steps, x_0 and y_0 is the initial value to start from, range is how long we want to go
		int computation_steps = (range / delta_x)+1;
		std::cout << computation_steps * 2 << std::endl;

		std::vector<float> values;
		Timer t1;
		t1.Begin();
		values.resize(2*computation_steps);
		std::cout << t1.getTimeNano() << std::endl;
		values.push_back(x_0);
		values.push_back(y_0);

		float last_y = y_0;
		float push = .1;
		int index = 0;
		for (int i = 0; i < computation_steps; i++)
		{
			float cc = y_derivative(i, last_y);
			values[index++] = cc;
			values[index++] = cc;
		}
		return values;
		for (float x = x_0+delta_x; x < x_0 + range; x += delta_x)
		{
			//t1.Begin();
			float new_y = last_y + y_derivative(x,last_y) * delta_x;
			//if (std::isinf(new_y))
			//{
				//were at infinity so take the last noninfinity and like reflect that
				//x = values[values.size() - 2] + push;
				//last_y will already be that
				//continue;
				//break;
			//}
			values.push_back(x);
			values.push_back(new_y);

			//std::cout << x << " " << new_y << std::endl;

			last_y = new_y;
			//std::cout<<t1.getTimeNano()<<std::endl;
		}

		std::vector<float> backwardvalues;
		backwardvalues.reserve(4 * computation_steps);
		last_y = y_0;
		for (float x = x_0-delta_x; x > x_0 - range; x -= delta_x)
		{
			float new_y = last_y - y_derivative(x, last_y) * delta_x;
			//if (std::isinf(new_y))
			//{
				//std::cout << "inf\n";
				//were at infinity so take the last noninfinity and like reflect that
				//x = values[values.size() - 1] - push;
				//last_y will already be that
				//continue;
				//break;
			//}
			backwardvalues.push_back(new_y);
			backwardvalues.push_back(x);

			//std::cout << x << " " << new_y << std::endl;

			last_y = new_y;
		}
		
		std::reverse(backwardvalues.begin(), backwardvalues.end());

		for (int i = 0; i < values.size(); i++)
		{
			backwardvalues.push_back(values[i]);
		}
		


		return backwardvalues;
	}

	//backwards euler
	/*
		so this method
	
	*/

	void MainLoop()
	{
		//initial conditions, boundary conditions
		//y(x0)=y0, or y'(x0)=y'0 or a boundary y(x0)=y0 and y(x1)=y1


		//Basic euler
		//y_n+1 = y_n * f(x,y)*x_step

		//RK2

		//RK4

		//multistep

		//also I can visualize each of the steps, and all the stuff now with my plotter and the mouse_to_world coordinate math

		glm::vec2 graph_scale(1850, 1850 / 1.5f);
		glm::vec2 graph_scale2(600.0f, 600.0 / 1.5f);


		//graph1.createsin(100, .05);

		//Graph2D graph2(graphics_manager, graph_scale2.x, graph_scale2.y);
		//graph2.setbackgroundcolor(glm::vec3(1, 1, 1));
		//graph2.createsin(100, .05);


		std::vector<const char*> functionmapname = {"y'=6+3y"};
		std::vector<const char*> methodmapname = { "Eulers method", "f(x)+f(x+2h)/2 - f'(x+h)-f'(x)/h" };

		int current_method = 0;
		int current_function = 0;

		float delta_x = .01f;
		float range = 100.f;
		float iv[2] = { 0,4 };
		bool show_1 = true;
		bool show_2 = false;
		bool enable_mouse_iv = false;

		Graph2D graph1(graphics_manager, graph_scale.x, graph_scale.y);
		graph1.setbackgroundcolor(glm::vec3(1, 1, 1));
		graph1.enableline(1, show_2);

		while (!glfwWindowShouldClose(graphics_manager.getWindow()))
		{
			//INPUT
			input_manager.HandleInput();
			glm::vec2 mouse_pos = graphics_manager.WindowToWorld(input_manager.getMousePosition());

			bool holding_not_first = (input_manager.getMouseButtonCurrent().x && !input_manager.getMouseButtonOnce().x);
			graph1.give_mouse_data(input_manager.getMousePosition(), input_manager.MouseScrollDirection(), holding_not_first, input_manager.getMouseButtonOnce().x);
			//graph2.give_mouse_data(input_manager.getMousePosition(), input_manager.MouseScrollDirection(), holding_not_first, input_manager.getMouseButtonOnce().x);

			glm::vec3 graph_translate(0, 0, 0);
			//glm::vec2 graph_scale(graphics_manager.getWidth(), graphics_manager.getHeight()); //x, x/1.5f
			graph1.setviewport(graph_translate.x, graph_translate.y, graph_scale.x, graph_scale.y);

			glm::vec3 graph_translate2(250, -250, 0);
			//graph2.setviewport(graph_translate2.x, graph_translate2.y, graph_scale2.x, graph_scale2.y);


			//graph predraw
			graph1.DrawToTexture();
			//graph2.DrawToTexture();

			//DRAW
			graphics_manager.Begin();

			//draw the plotter texture
			glm::mat4 mat = glm::translate(glm::mat4(1.0f), graph_translate);
			mat = glm::scale(mat, glm::vec3(graph_scale.x, graph_scale.y, 1)); //3:2 ratio
			DrawSubmit submit_1{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::TEXTURE, mat, glm::vec3(50,0,0) };
			submit_1.sampler_id = graph1.gettexture();
			graphics_manager.SubmitDraw(submit_1);

			//glm::mat4 mat2 = glm::translate(glm::mat4(1.0f), graph_translate2);
			//mat2 = glm::scale(mat2, glm::vec3(graph_scale2.x, graph_scale2.y, 1)); //3:2 ratio
			//DrawSubmit submit_2{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::TEXTURE, mat2, glm::vec3(50,0,0) };
			//submit_2.sampler_id = graph2.gettexture();
			//graphics_manager.SubmitDraw(submit_2);

			graphics_manager.Render();

			glm::vec2 zzz = graph1.mousetographpos(input_manager.getMousePosition());
			std::string zzzz = std::to_string(zzz.x) + ", " + std::to_string(zzz.y);
			graphics_manager.RenderTextFont(zzzz, 450, 450, 0.5f, glm::vec3(0, 0, 0));

			//IMGUI
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();

			ImGui::NewFrame();
			//ImGui::ShowDemoWindow();
			ImGui::Begin("ODE Solver app");

			ImGui::Combo("solver methods", &current_method, methodmapname.data(), methodmapname.size(), 15);
			ImGui::Combo("functions", &current_function, functionmapname.data(), functionmapname.size(), 5);

			ImGui::InputFloat("delta_x", &delta_x,.1,.1,4);
			ImGui::SliderFloat("range", &range, 1, 10000);
			ImGui::InputFloat2("initial value (x,y)", iv);

			if (ImGui::Button("Calculate") || enable_mouse_iv)
			{
				//y=-2+ce^3x
				auto derivative_1 = [=](float x, float y) -> float { return 6 + 3 * y; };
				auto solution_1 = [=](float x) -> float { return -2+((iv[1] + 2) / std::exp(3 * iv[0])) * std::exp(3 * x); };

				//y=x-1+ce^-x
				auto derivative_2 = [=](float x, float y) -> float { return x-y; };
				auto solution_2 = [=](float x) -> float { return x - 1 + ((iv[1] - iv[0] + 1) / std::exp(-iv[0])) * std::exp(-x); };
				
				//y^2(1-x^2)-cos^2x=c
				auto derivative_3 = [=](float x, float y) -> float { return (x*y*y - std::cos(x)*std::sin(x)) / (y*(1-x*x)); };
				//auto solution_3 = [=](float x) -> float { return x - 1 + ((iv[1] - iv[0] + 1) / std::exp(-iv[0])) * std::exp(-x); };

				auto derivative_4 = [=](float x, float y) -> float { return (x * (1 - x)) / (y*(-2+y)); };
				//auto solution_3 = [=](float x) -> float { return x - 1 + ((iv[1] - iv[0] + 1) / std::exp(-iv[0])) * std::exp(-x); };

				auto derivative_5 = [=](float x, float y) -> float { return std::sin(y); };
				
				auto derivative_6 = [=](float x, float y) -> float { return std::sqrt(x*y); };
				
				std::array<std::vector<float>,10> points;
				
				glm::vec2 theiv = glm::vec2(iv[0], iv[1]);
				if (enable_mouse_iv)
					theiv = zzz;
				switch (current_method)
				{
				case 0:
					//eulers
					for (int i = 0; i < 10; i++)
					{
						if (i == 0)
						{
							points[i] = eulersmethod(delta_x, theiv.x, theiv.y, range, derivative_1);
						}
						else
						{
							//points[i] = eulersmethod(delta_x, -50 + i*5, -50 + 100*std::pow(-1,i), range, derivative_5);
						}
					}
					break;
				}

				graph1.clearline(0);
				graph1.setlinecolor(glm::vec3(255, 0, 0), 0);
				graph1.plotpoints(0, points[0]);
				
				graph1.clearline(1);
				graph1.setlinecolor(glm::vec3(200, 200, 0), 1);
				graph1.createfunction(1, 1000, .1, solution_2);

				for (int i = 2; i < 10; i++)
				{
					//graph1.clearline(i);
					//graph1.setlinecolor(glm::vec3(255, 0, 0), 0);
					//graph1.plotpoints(i, points[i]);
				}
			}

			if (ImGui::Checkbox("show line 1 (euler)", &show_1))
			{
				graph1.enableline(0, show_1);
			}
			if (ImGui::Checkbox("show line 2 (real solution)", &show_2))
			{
				graph1.enableline(1, show_2);
			}
			if (ImGui::Checkbox("make mouse initial point", &enable_mouse_iv))
			{
			
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