#pragma once
#include "../Graphics/GraphicsManager.h"
#include "../Graphics/Graph2D.h"

class PlotterTestApp
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		PlotterTestApp* handler = reinterpret_cast<PlotterTestApp*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<PlotterTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<PlotterTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<PlotterTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<PlotterTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons;
		if (graphics_manager.InitWindow(1850, 1850/1.5f, "PlotterTest App", icons) == EXIT_FAILURE)
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
		//Basic euler
		//y_n+1 = y_n * f(x,y)*x_step

		//RK2

		//RK4

		//multistep

		//also I can visualize each of the steps, and all the stuff now with my plotter and the mouse_to_world coordinate math

		glm::vec2 graph_scale(1850, 1850 / 1.5f);
		glm::vec2 graph_scale2(600.0f, 600.0 / 1.5f);

		Graph2D graph1(graphics_manager, graph_scale.x, graph_scale.y);
		graph1.setbackgroundcolor(glm::vec3(1, 1, 1));
		graph1.createsin(100, .05);

		Graph2D graph2(graphics_manager, graph_scale2.x, graph_scale2.y);
		graph2.setbackgroundcolor(glm::vec3(1, 1, 1));
		graph2.createsin(100, .05);

		while (!glfwWindowShouldClose(graphics_manager.getWindow()))
		{
			//INPUT
			input_manager.HandleInput();
			glm::vec2 mouse_pos = graphics_manager.WindowToWorld(input_manager.getMousePosition());
			
			bool holding_not_first = (input_manager.getMouseButtonCurrent().x && !input_manager.getMouseButtonOnce().x);
			graph1.give_mouse_data(input_manager.getMousePosition(), input_manager.MouseScrollDirection(), holding_not_first, input_manager.getMouseButtonOnce().x);
			graph2.give_mouse_data(input_manager.getMousePosition(), input_manager.MouseScrollDirection(), holding_not_first, input_manager.getMouseButtonOnce().x);

			glm::vec3 graph_translate(0, 0, 0);
			//glm::vec2 graph_scale(graphics_manager.getWidth(), graphics_manager.getHeight()); //x, x/1.5f
			graph1.setviewport(graph_translate.x, graph_translate.y, graph_scale.x, graph_scale.y);
			
			glm::vec3 graph_translate2(250, -250, 0);
			graph2.setviewport(graph_translate2.x, graph_translate2.y, graph_scale2.x, graph_scale2.y);
			
			
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

			glm::mat4 mat2 = glm::translate(glm::mat4(1.0f), graph_translate2);
			mat2 = glm::scale(mat2, glm::vec3(graph_scale2.x, graph_scale2.y, 1)); //3:2 ratio
			DrawSubmit submit_2{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::TEXTURE, mat2, glm::vec3(50,0,0) };
			submit_2.sampler_id = graph2.gettexture();
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
			ImGui::Begin("test app");
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