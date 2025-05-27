#pragma once
#include "../Graphics/GraphicsManager.h"
#include <cerrno>
#include <cfenv>

class MuscleCrossbridgeApp
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		MuscleCrossbridgeApp* handler = reinterpret_cast<MuscleCrossbridgeApp*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<MuscleCrossbridgeApp*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<MuscleCrossbridgeApp*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<MuscleCrossbridgeApp*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<MuscleCrossbridgeApp*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons = {
			"Game_Data/muscle16.png",
			"Game_Data/muscle32.png",
			"Game_Data/muscle48.png"
		};
		if (graphics_manager.InitWindow(1500, 1000, "Muscle Crossbridge Simulation", icons) == EXIT_FAILURE)
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

	//v(x) functions
	double velocity_spring(double x)
	{
		int direction = cos(total_time);

		//return A * sin(1 * x);
		float pos = std::clamp(x, -A, A);
		float eps = 0.0f;
		if (abs(pos-A) <= .00001)
			eps = .000001;
		float k = v;// 150;
		
		//std::cout << k * A * std::sin(std::acos((x- eps) / A)) << std::endl;
		
		if (pos > 0)
			return k * A * std::sin(std::acos((x - eps) / A));
		else
		{
			if (abs(pos + A) <= .5)
				return 0;
			return k * A * std::sin(std::acos((x) / A));
		}
	}

	double velocity_global_oscillator(double t)
	{
		float k = 5.0;
		return v * sin(2 * 3.14 * k * t);
	}

	//f(x) functions
	//force of a crossbridge given the distance from equilibrium (can be positive or negative
	double force_of_bridge(float x)
	{
		errno = 0;
		std::feclearexcept(FE_ALL_EXCEPT);
		double exponent = std::exp(gamma * x);
		if (errno == ERANGE || std::fetestexcept(FE_OVERFLOW))
		{
			std::cout << "exponent too big\n";
			return 0;
		}

		return p1*(exponent - 1);
	}
	double force_of_bridge_linear(float x)
	{
		return p1*x;
	}

	//VARIABLES
	/*
	length = nanometers (10^-9)
	time = second
	velocity = nm/s
	force = pN (picoNewtons 10^-12 N)
	*/
	uint32_t n_0 = 10000; //number of crossbridges (unitless)
	uint32_t N = 6000; //number of sarcomeres (unitless)

	double atcht_rate = 14; //attatchment rate (probability / unit time)
	double det_rate = 126; //dettatchment rate (probability / unit time)

	std::vector<uint32_t> attatch_info;
	std::vector<double> crossbridgesposition;

	double total_time = 0.0; // (seconds)
	double Time_step = 0.01 / (atcht_rate + det_rate); // (seconds)

	double A = 5; //starting attatchment distance (length)
	double gamma = 0.322; //(distance)
	double p1 = 4; //(force)

	double P = 0; //force (mass*length/time^2)
	double v = 500; //velocity of half sarcamere (distance / time)
	double V = v * (2 * N); //total velocity of muscle (length / time)

	int velocity_mode = 0;
	float still_time = total_time;
	int time_count = 0;

	void MainLoop()
	{
		//This is app is to simulate and explore the muscle crossbridge theory related to the force-velocity curve
		//Resource: Modeling and Simulation in Medicine and the Life Sciences, Fran C. Hoppensteadt and Charles S.Peskin, Chapter 5

		//U(X)
		//first we need the population density function u(x), every attatched crossbridge has some distance x from equilibrium, u(x) gives
		//the percentage of total crossbridges that are at a distance x. Close to a probability distribution function but some can be detatched so it doesnt
		//sum to 1

		//instead of analytically computing this, we can just simulate it. We attatch, detatch crossbridges and move them with velocity and
		//u(x) is simply creating a type of histogram where every crossbridges current x value goes into a specified bin, we can compute the discrete data set

		//P
		//next P is the force of half of 1 sarcomere where the crossbridges live, this relates to the force the muscle feels. Given a p(x), this value is simply
		//going over every crossbridge and gettings its force p(x), then adding all those up, this will create the actual force at a given time

		//time
		//so the only difference in the theory is that we have time now. In our simulation we are simulating everything with time steps, and u(x) also its just 
		//steady state (which we assume), but it will eventually get steady state after some time. 
		
		//random generation
		//also we will be using random numbers for the attatchement and detratchement which the theory assumes are constant rates based off of how many are in the
		//detatched or attatched pool, it will just be uniform random number times the rate A or B of attach/detatch
		//We say the probability with time rate is v, so if some timestep happens generally the probability is v*timestep, this is an approximation of course, and
		//as the timestep gets bigger the approximation gets worse and worse, its best if timestep is smallest and approaching 0. Though in our theory a and b rates
		//are constant so it should not really matter, the change of probability for a change of time is just constant based off the slope

		//Problem 1: given the whole muscles shortening time as a function of time V(t), what is the corresponding force P?

		//Problem 2: given the whole muscles force with  time P(t), what is the corresponding motion and speed V(t)?

		if (abs(force_of_bridge(0) / force_of_bridge(A)) > 10 * std::numeric_limits<float>::min())
			std::cout << "precision error!\n";

		bool paused = false;

		attatch_info.clear();
		crossbridgesposition.clear();
		attatch_info.resize(n_0, 0);
		crossbridgesposition.resize(n_0, 0);

		//std::vector<std::pair<float, float>> velocity_force_points;
		std::vector<float> velocity_force_points;

		int n_0_tmp = n_0;
		int N_tmp = N;
		float a_tmp = atcht_rate;
		float b_tmp = det_rate;
		float t_tmp = Time_step;
		float A_tmp = A;
		float v_tmp = v;
		float p1_tmp = p1;
		float g_tmp = gamma;
		int vel_mod_tmp = velocity_mode;

		if ((atcht_rate * Time_step) > 1 || (det_rate * Time_step) > 1)
			std::cout << "probability greater than 1\n";

		bool imgui_hovered = false;
		while (!glfwWindowShouldClose(graphics_manager.getWindow()))
		{
			//INPUT
			input_manager.HandleInput();
			glm::vec2 mouse_pos = graphics_manager.WindowToWorld(input_manager.getMousePosition());

			//SIMULATION
			if (!paused)
			{
				Simulate(false, 0);
			}
			
			
			//DRAW
			graphics_manager.Begin();

			graphics_manager.submitPlot(0);

			graphics_manager.Render();

			//IMGUI
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();
			
			ImGui::NewFrame();
			//ImGui::ShowDemoWindow();
			ImGui::Begin("Muscle CrossBridge Simulation");	
			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

			//VARIABLES
			ImGui::SliderInt("(n0) crossbridge count", &n_0_tmp, 0, 30000);
			ImGui::SliderInt("(N) sarcamere count in whole muscle", &N_tmp, 0, 10000);
			ImGui::SliderFloat("(a) attach rate", &a_tmp, 0, 20);
			ImGui::SliderFloat("(b) detach rate", &b_tmp, 0, 150);
			ImGui::SliderFloat("(A) starting distance", &A_tmp, 0, 5);
			ImGui::SliderFloat("(v) sarcamere muscle velocity", &v_tmp, 0, 2000);
			ImGui::SliderFloat("(p1) p(x) constant", &p1_tmp, 0, 5);
			ImGui::SliderFloat("(gamma) p(x) constant", &g_tmp, 0, 1);
			ImGui::SliderFloat("(dt) timestep", &t_tmp, 0,.001,"%.8f");

			std::string s1 = "constant velocity";
			std::string s2 = "spring velocity";
			std::string s3 = "global oscillator";
			std::array<const char*, 3> modes = { s1.c_str(), s2.c_str(), s3.c_str()};
			ImGui::Combo("velocity mode", &vel_mod_tmp, modes.data(),3);

			std::string val_10 = "Total Simulation Time Elapsed: " + std::to_string(total_time);
			ImGui::Text(val_10.c_str());

			if (time_count % 100 == 0)
			{
				still_time = P;
			}
			time_count = (time_count + 1) % 100000;

			std::string val_11 = "force (P): " + std::to_string(still_time);
			ImGui::Text(val_11.c_str());

			std::string val_12 = "v: " + std::to_string(v);

			ImGui::ProgressBar(P / 20000, ImVec2(0, 15), "P");
			ImGui::ProgressBar(v/2000, ImVec2(0, 15), val_12.c_str());

			if (ImGui::Button("Restart Simulation"))
			{
				n_0 = n_0_tmp;
				N = N_tmp;
				atcht_rate = a_tmp;
				det_rate = b_tmp;
				A = A_tmp;
				Time_step = t_tmp;
				total_time = 0.0;
				v = v_tmp;
				V = v * (2 * N);
				p1 = p1_tmp;
				gamma = g_tmp;
				P = 0.0;
				velocity_mode = vel_mod_tmp;


				if (abs(force_of_bridge(0) / force_of_bridge(A)) > 10 * std::numeric_limits<float>::min())
					std::cout << "precision error!\n";

				attatch_info.clear();
				crossbridgesposition.clear();
				attatch_info.resize(n_0, 0);
				crossbridgesposition.resize(n_0, 0);
			}
			if (!paused)
			{
				if (ImGui::Button("Pause Simulation"))
				{
					paused = true;
				}
			}
			else
			{
				if (ImGui::Button("Play Simulation"))
				{
					paused = false;
				}
			}

			//HISTOGRAM DATA
			const uint32_t bin_count = 20;
			const float beggining_x = -2*A;
			const float end_x = A;
			const float delta_x = (end_x - beggining_x) / (float)bin_count;

			std::array<float, bin_count> histogram;
			std::array<float, bin_count> histogram_ratios;
			histogram.fill(0);
			histogram_ratios.fill(0);
			for (int i = 0; i < n_0; i++)
			{
				if (attatch_info[i] == 1)
				{
					if (crossbridgesposition[i] < beggining_x) {
						//std::cout << "position less than -A "<< crossbridgesposition[i]<< "\n";
						continue;
					}
				
					uint32_t bin_number = std::floor(std::clamp((crossbridgesposition[i] - beggining_x) / delta_x,0.0,(double)(bin_count-1)));
					histogram[bin_number] += 1;
					histogram_ratios[bin_number] += 1;
				}
			}
			for (int i = 0; i < histogram_ratios.size(); i++)
			{
				//just divide by total number so every value is a ratio
				histogram_ratios[i] /= (float)n_0;
			}

			std::array<float, 2> attatch_plot;
			attatch_plot.fill(0);
			for (int i = 0; i < attatch_info.size(); i++)
			{
				bool val = attatch_info[i];
				attatch_plot[val]++;
			}

			std::string val = "u(x) [" + std::to_string(beggining_x) + "," + std::to_string(end_x) + "]";
			std::string val2 = "# of crossbridges in interval\nsize: " + std::to_string(delta_x);
			//ImGui::PlotHistogram(val2.c_str(), histogram.data(), histogram.size(), 0, val.c_str(), 0.0f, .2*n_0, ImVec2(0, 180));

			std::string val6 = "% of crossbridges in interval\nsize: " + std::to_string(delta_x);
			ImGui::PlotHistogram(val6.c_str(), histogram_ratios.data(), histogram_ratios.size(), 0, val.c_str(), 0.0f, 0.2f, ImVec2(0, 180));
			
			ImGui::PlotHistogram("", attatch_plot.data(), attatch_plot.size(), 0, "U(x) detatched vs attatched", 0.0f, n_0, ImVec2(0, 100));

			std::vector<float> float_vers(n_0);
			for (int i = 0; i < crossbridgesposition.size(); i++)
			{
				float_vers[i] = static_cast<float>(crossbridgesposition[i]);
			}

			std::string val3 = std::to_string(end_x) + "\n\n\n\n\n\n\n\n\n\n\n\n\n" + std::to_string(beggining_x);
			ImGui::PlotHistogram(val3.c_str(), float_vers.data(), float_vers.size(), 0, "position data", -A, A, ImVec2(0, 180));


			if (ImGui::Button("Calculate Force-Velocity Curve"))
			{
				//pick different values of v, and get the equilibrium P. Then we know what P is and V=v2N
				//for v[0,1000] what is the corresponding P? over .1 seconds equilibrium time

				//run simulation with v value x for .1 seconds, then give the P value
				float original_v = v;
				velocity_force_points.clear();
				velocity_force_points.resize(2300);
				int precision = 10;
				for (int i = 0; i < velocity_force_points.size(); i++)
				{
					if (i % precision == 0)
					{
						v = i; //set velocity
						float P_val = Simulate(true, .01);
						std::cout << P_val << std::endl;
						velocity_force_points[i] = P_val;
					}
					else
					{
						velocity_force_points[i] = velocity_force_points[i - (i % precision)];
					}
				}
				//std::reverse(velocity_force_points.begin(), velocity_force_points.end());
				v = original_v;

				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0,0,0));
				mat = glm::scale(mat, glm::vec3(.05, .1, 1));

				graphics_manager.CreatePlot(velocity_force_points.data(), velocity_force_points.size(), mat);
			}
			std::string val11 = std::to_string(20000) + "\n\n\n\n\n\n\n\n\n\n\n\n\n" + std::to_string(0);
			ImGui::PlotLines(val11.c_str(), velocity_force_points.data(), velocity_force_points.size(), 0, "Force-Velocity curve", 0, 20000, ImVec2(0, 180));

			//ImGui::begin

			ImGui::End();
			
			ImGui::Render();
			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

			graphics_manager.End();
		}
	}

	double Simulate(bool start, float time_to_simulate)
	{
		if (start)
		{
			attatch_info.clear();
			crossbridgesposition.clear();
			attatch_info.resize(n_0, 0);
			crossbridgesposition.resize(n_0, 0);
			total_time = 0;
		}

		do
		{
			total_time += Time_step;

			double a_probability = std::clamp(atcht_rate * Time_step, 0.0, 1.0);
			double d_probability = std::clamp(det_rate * Time_step, 0.0, 1.0);
			for (int i = 0; i < n_0; i++)
			{
				if (attatch_info[i] == 0)
				{
					double num = Helper::GetRandomNumberD(0.0, 1.0);
					if (a_probability >= num)
					{
						attatch_info[i] = 1;
						crossbridgesposition[i] = A;
					}
				}
				else
				{
					double num = Helper::GetRandomNumberD(0.0, 1.0);
					if (d_probability >= num) {
						attatch_info[i] = 0;
						crossbridgesposition[i] = 0;
					}
					else
					{
						//for constant
						if(velocity_mode == 0)
							crossbridgesposition[i] -= v * Time_step;
						
						//for spring
						if (velocity_mode == 1)
						{
							float vel = velocity_spring(crossbridgesposition[i]);
							if (vel == 0)
							{
								//were at the max stretching on -A so just break it
								attatch_info[i] = 0;
								crossbridgesposition[i] = 0;
							}
							else
							{
								crossbridgesposition[i] -= vel * Time_step;
							}
						}

						if (velocity_mode == 2)
						{
							crossbridgesposition[i] -= velocity_global_oscillator(total_time) * Time_step;
						}
					}
				}
			}

			//for the total force, just add up every crossbridges force component, which is just p(x)
			P = 0;
			for (int i = 0; i < n_0; i++)
			{
				if (attatch_info[i] == 1)
					P += force_of_bridge(crossbridgesposition[i]);
			}

			//std::cout << crossbridgesposition[3] << std::endl;
		} while (total_time <= time_to_simulate);

		return P;
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