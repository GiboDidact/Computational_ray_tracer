#pragma once
#include "../../Graphics/GraphicsManager.h"
#include "../../RayTracer/MonteCarlos.h"
class MonteCarlosTestApp
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		MonteCarlosTestApp* handler = reinterpret_cast<MonteCarlosTestApp*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<MonteCarlosTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<MonteCarlosTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<MonteCarlosTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<MonteCarlosTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons;
		if (graphics_manager.InitWindow(1850, 1850 / 1.5f, "MonteCarlos App", icons) == EXIT_FAILURE)
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
		//write basic monte carlos integrator
		//look at the convergence rate and find variance
		//do another sampling distribution 
		
		//write biased one
		//do a variance reduction technique
		//do some sampling on a data set
		//multivariable integrator

		float a = 5;
		float b = 12;

		auto function_1 = [=](float x) -> float { return 6.0f + 3.0f*x; };
		float function_1_answer = 6 * b - 6 * a + (3.0f / 2.0f) * (b * b - a * a);

		auto function_2 = [=](float x) -> float { return cos(x)+5; };
		float function_2_answer = 35.4223513567; //(5,12)

		float U = Helper::GetRandomNumber(0.0, 1.0); //~uniform(0,1)
		

		float n = 1;
		float sum = 0;
		float estimator_1_value = 0.0f;
		float sum_2 = 0;
		float estimator_2_value = 0.0f;

		float MSE_total = 0.0f;
		float MSE = 0.0f;

		auto test_function = function_2;
		float test_answer = function_2_answer;

		float delta = 0.1f;
		auto exponential_dist = [=](float x) -> float { return delta * std::expf(-delta * (x - a)); };
		auto exponential_dist_normalized = normalizepdf(exponential_dist, a, b);
		auto exponential_dist2 = [=](float x) -> float { return (delta) * std::expf(delta * (x - a)); };
		auto exponential_dist2_normalized = normalizepdf(exponential_dist2, a, b);

		auto linear_dist = [=](float x)-> float { return 6 + 3 * x; };
		auto lienar_dist_normalized = normalizepdf(linear_dist, a, b);

		float zz = 3.5f;
		auto normal_dist = [=](float x) -> float { return (1.0f / (std::abs(zz) * sqrt(std::numbers::pi))) * std::exp(-std::powf((x-8.5f) / zz, 2.0f)); };
		auto normal_dist_normalized = normalizepdf(normal_dist, a, b);

		uniformEstimator estimator1;
		distributionEstimator estimator2;
		float EV = estimator2.approximate_expected_value(test_function, normal_dist_normalized, a, b);
		float VAR = estimator2.approximate_variance(test_function, normal_dist_normalized, a, b);
		float EFF = estimator2.approximate_efficiency(test_function, normal_dist_normalized, a, b);
		estimator1.PrintChebychevInequalityRanges(test_function, a, b);

		estimator_1_value = estimator2.EstimateFunctionN(100, test_function, normal_dist_normalized, a, b);
		estimator_2_value = estimator2.EstimateFunctionN(10000, test_function, normal_dist_normalized, a, b);

		//test the max estimator of b, we have a uniform distribution of (b-a), we want to estimate b with max(x1,x2...,xn)
		//the longterm bias is a, so the expected error in the long term is a, else is n/n+1a -b/n+1
		/*int N = 5000;
		float themax = std::numeric_limits<float>::min();
		for (int i = 0; i < N; i++)
		{
			float X_i = Helper::GetRandomNumber(a, b);
			themax = glm::max(themax, X_i);
		}
		std::cout <<N<<" samples estimation of b: " << b << " is: " << themax << std::endl;
		*/
		while (!glfwWindowShouldClose(graphics_manager.getWindow()))
		{
			//INPUT
			input_manager.HandleInput();
			glm::vec2 mouse_pos = graphics_manager.WindowToWorld(input_manager.getMousePosition());


			//DRAW
			graphics_manager.Begin();
			graphics_manager.Render();

			//IMGUI
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();

			ImGui::NewFrame();
			//ImGui::ShowDemoWindow();
			ImGui::Begin("monte carlos app");

			std::vector<float> data;
			data.push_back(test_answer);
			data.push_back(estimator_1_value);
			data.push_back(estimator_2_value);
			ImGui::PlotHistogram("integrator", data.data(), data.size(), 0, 0, 0.0f, 350.0f, ImVec2(0, 180));
			
			std::string s1 = "n: " + std::to_string(n);
			ImGui::Text(s1.c_str());

			std::string s2 = "MSE: " + std::to_string(MSE);
			ImGui::Text(s2.c_str());

			std::string s3 = "Expected Value: " + std::to_string(EV);
			ImGui::Text(s3.c_str());

			std::string s4 = "Variance: " + std::to_string(VAR);
			ImGui::Text(s4.c_str());

			std::string s5 = "Efficiency: " + std::to_string(EFF);
			ImGui::Text(s5.c_str());

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