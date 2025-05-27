#pragma once
#include "../../Graphics/GraphicsManager.h"
#include "../../Graphics/Graph2D.h"
#include "../../RayTracer/Sampling.h"
#include "../../ThirdParty/pbrv4/samplers.h"
using namespace pbrt;

class SamplerTestApp
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		SamplerTestApp* handler = reinterpret_cast<SamplerTestApp*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<SamplerTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<SamplerTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<SamplerTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<SamplerTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons;
		if (graphics_manager.InitWindow(1850, 1850 / 1.5f, "Sampler Test App", icons) == EXIT_FAILURE)
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
		std::vector<const char*> sampling_methods = { "continuous_inversion", "continuous_rejection", "discrete_inversion", "discrete_rejection",
		"analytic_linear", "analytic_tent", "analytic_exponential", "analytic_normal", "analytic_logistic", "analytic_bilinear",
		"uniform_disk_polar", "uniform_disk_concentric", "uniform_hemisphere", "uniform_sphere", "cosine_hemisphere", "cone",
		"uniform pixel samplers"};
		int current_method = 16;
		int previous_method = -1;
		const char* name = "sampling methods";

		//Generate all test case data
		
		//CONTINUOUS (inversion/rejection)
		float delta = 2.0f;
		auto exponential_dist = [=](float x) -> float { return delta * std::expf(-delta * x); };
		Continuous_Inversion_Sampler exponentialsampler(exponential_dist, 0, 2, 1000);
		std::vector<float>exp_samples = exponentialsampler.generateTestSampleHistogram(100, 1000000, 25.0f, -1);
		std::vector<float>exp_sample_rejection = generateContinuousRejectionTestSampleHistogram(exponential_dist, 0, 3, 100, 10000, 25.0f, -1);
		
		auto xx_dist = [=](float x) -> float { return (x * x) / 9.0; };
		Continuous_Inversion_Sampler xsquaredsampler(xx_dist, 0, 3, 50);
		Timer t1;
		std::vector<float>xx_samples = xsquaredsampler.generateTestSampleHistogram(100, 100000, 25.0f, -1);
		std::vector<float>xx_rejection_histogram = generateContinuousRejectionTestSampleHistogram(xx_dist, 0, 3, 100, 100000, 25.0f, -1);
	
		float zz = .2f;
		auto normal_dist = [=](float x) -> float { return (1.0f/(std::abs(zz)*sqrt(std::numbers::pi)))*std::exp(-std::powf(x / zz,2.0f)); };
		Continuous_Inversion_Sampler normal_sampler(normal_dist, -2, 2, 1000);
		std::vector<float>normal_samples = normal_sampler.generateTestSampleHistogram(100, 100000, 25.0f, 0);
		std::vector<float>normal_rejection_histogram = generateContinuousRejectionTestSampleHistogram(normal_dist, -2, 2, 100, 100000, 25.0f, -1);

		//DISCRETE (inversion/rejection)
		std::vector<float> geometric_dist;
		std::vector<float> geom_xvalues;
		for (int i = 1; i < 25; i++)
		{
			float p = .05;
			geometric_dist.push_back(std::powf(1 - p, i - 1) * p);
			geom_xvalues.push_back(i);
		}
		Discrete_Inversion_Sampler geometric_sampler(geometric_dist);
		std::vector<float> geometric_histogram = geometric_sampler.generateTestSampleHistogram(geom_xvalues, 100000, 50.0f);
		std::vector<float> geometric_rejection_histogram = generateRejectionTestSampleHistogram(geom_xvalues, geometric_dist, 100000, 50.0f);

		std::vector<float> binomial_dist;
		std::vector<float> binomial_xvalues;
		float trials = 25;
		for (int successes = 0; successes < trials; successes++)
		{
			float p = .2;
			binomial_dist.push_back(Helper::nCr(trials, successes) * std::powf(p, successes) * std::powf(1 - p, trials - successes));
			binomial_xvalues.push_back(successes);
		}
		Discrete_Inversion_Sampler binomial_sampler(binomial_dist);
		std::vector<float> binomial_histogram = binomial_sampler.generateTestSampleHistogram(binomial_xvalues, 10000, 10);
		std::vector<float> binomial_rejection_histogram = generateRejectionTestSampleHistogram(binomial_xvalues, binomial_dist, 10000, 10);
		
		//sample linear, (1-x)a + xb, a,b>=0, x[0,1]y[a,b]
		float lin_a = 5;
		float lin_b = 1;
		auto linear_func = [=](float x) -> float { return (1 - x) * lin_a + x * lin_b; };
		auto linear_func_normalized = normalizepdf(linear_func, 0, 1);
		std::vector<float> linear_samples;
		linear_samples.reserve(1000000);
		for (int i = 0; i < 1000000; i++)
		{
			float X_i = SampleLinear(lin_a, lin_b);
			linear_samples.push_back(X_i);
		}
		std::vector<float> linear_histogram = generateContinuousSampleHistogram(linear_samples, 0, 1, 100, 100.0f, -1);

		//sample tent
		float tent_r = 2.0f;
		auto tent_func = [=](float x) -> float { if (std::abs(x) < tent_r) { return tent_r - std::abs(x); } return 0; };
		auto tent_func_normalized = normalizepdf(tent_func, -tent_r, tent_r);
		std::vector<float> tent_samples;
		tent_samples.reserve(10000);
		for (int i = 0; i < 100000; i++)
		{
			float X_i = SampleTent(tent_r);
			tent_samples.push_back(X_i);
		}
		std::vector<float> tent_histogram = generateContinuousSampleHistogram(tent_samples, -tent_r, tent_r, 100, 100.0f, -1);

		//sample exponential
		float exp_a = 0.5f; //can't be 0
		float x_range_1 = (exp_a < 0) ? -10.0f : 0;
		float x_range_2 = (exp_a > 0) ? 10.0f : 0;
		auto exp_func_normalized = [=](float x) -> float { return std::abs(exp_a)*std::exp(-exp_a*x); };
		std::vector<float> exp2_samples;
		exp2_samples.reserve(10000);
		for (int i = 0; i < 100000; i++)
		{
			float X_i = SampleExponential(exp_a);
			exp2_samples.push_back(X_i);
		}
		std::vector<float> exp_histogram = generateContinuousSampleHistogram(exp2_samples, x_range_1, x_range_2, 100, 25.0f, -1);

		//gaussian
		float gauss_u = 0.0f;
		float gauss_s = 0.4f;
		float gauss_xrange_1 = gauss_u - 3 * gauss_s;
		float guass_xrange_2 = gauss_u + 3 * gauss_s;
		auto gauss_function_normalized = [=](float x) -> float { return (1.0f / std::sqrtf(2 * std::numbers::pi * gauss_s* gauss_s)) *
																		std::exp(-(std::pow(x - gauss_u, 2.0f)) / (2 * gauss_s* gauss_s)); };
		std::vector<float> gauss_samples;
		gauss_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			float X_i = SampleNormal(gauss_u, gauss_s);
			gauss_samples.push_back(X_i);
		}
		std::vector<float> gauss_histogram = generateContinuousSampleHistogram(gauss_samples, gauss_xrange_1, guass_xrange_2, 100, 25.0f, 0);

		//logistics
		float logistic_s = .3f;
		auto logistics_function_normalized = [=](float x) -> float { return std::exp(-std::abs(x) / logistic_s) / 
			                                                         (logistic_s *std::pow(1+ std::exp(-std::abs(x) / logistic_s), 2)); };
		std::vector<float> logistic_samples;
		logistic_samples.reserve(1000);
		for (int i = 0; i < 100000; i++)
		{
			float X_i = SampleLogistic(logistic_s);
			logistic_samples.push_back(X_i);
		}
		std::vector<float> logistic_histogram = generateContinuousSampleHistogram(logistic_samples, -5, 5, 100, 15.0f, 0);

		//bilinear, x[0,1], y[0,1]
		std::vector<float> bi_w = { 10,1,1,1 };
		auto bilinear_function = [=](float x, float y) {return (1 - x) * (1 - y) * bi_w[0] + x * (1 - y) * bi_w[1] + y * (1 - x) * bi_w[2] + x * y * bi_w[3]; };
		std::vector<glm::vec2> bilinear_samples;
		bilinear_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			glm::vec2 P = SampleBilinear(bi_w);
			bilinear_samples.push_back(P);
		}

		//uniform Disk
		std::vector<glm::vec2> diskpolar_samples;
		diskpolar_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			//SampleDiskNaive();
			//RejectionSampleDisk();
			//SampleUniformDiskPolar();
			glm::vec2 P = SampleUniformDiskPolar();
			diskpolar_samples.push_back(P);
		}
		IndependentSampler* indep_sampler2 = new IndependentSampler(1);
		indep_sampler2->StartPixelSample(glm::ivec2(0, 0), 0, 0);
		//uniform concentric Disk
		std::vector<glm::vec2> diskconcentric_samples;
		diskconcentric_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			glm::vec2 P = SampleUniformDiskConcentric(indep_sampler2->Get2D());
			diskconcentric_samples.push_back(P);
		}

		//uniform hemisphere
		std::vector<glm::vec2> hemisphere_samples;
		hemisphere_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			glm::vec2 P = SampleUniformHemisphere();
			hemisphere_samples.push_back(P);
		}

		//uniform sphere
		std::vector<glm::vec3> sphere_samples;
		sphere_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			glm::vec3 P = SampleUniformSphere();
			sphere_samples.push_back(P);
		}

		//cosine hemisphere
		std::vector<glm::vec3> cosinhemisphere_samples;
		cosinhemisphere_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			glm::vec3 P = SampleCosineHemisphere(indep_sampler2->Get2D());
			cosinhemisphere_samples.push_back(P);
		}

		//cone
		float costheta_max = std::cos(glm::radians(20.0f));
		std::vector<glm::vec3> cone_samples;
		cone_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			glm::vec3 P = SampleUniformCone(costheta_max);
			cone_samples.push_back(P);
		}

		/* SAMPLERS
			
			to test the samplers I want to have a grid of pixels and indices per pixel, then I do the sampling over all the pixels and indices
			the indices can have different colors 
			then to test determinism honestly I can jsut repeat a bunch of times and they should all be on top of eachother

			then to test dimensionality im not sure yet, I can figure that out later

		*/
		float sampler_scaling = 1000.0f;
		float sampler_point_scaling = 2.0f;
		std::vector<const char*> sampler_types = { "independant", "stratified", "stratified jitter", "sobol"};
		int sampler_type = 0; //0-independant, 1-stratified, 2-stratified_jitter
		bool sampler_everypixel = true;
		bool sampler_colorcoat = false;
		bool sampler_same_pixels = false;
		//x,y,indices counts

		const int X_pixels = 5;
		const int Y_pixels = 5;
		int pixel_indices = 10;
		std::array<std::array<std::vector<glm::vec2>, Y_pixels>, X_pixels>  indep_grid; //X by Y grid that holds an array of index vec2's

		//independant sampler
		int seed = 100;
		IndependentSampler* indep_sampler = new IndependentSampler(pixel_indices, seed);
		StratifiedSampler* strat_sampler = new StratifiedSampler(std::ceil(pixel_indices / 1.0f), std::ceil(pixel_indices / 1.0f), false, seed);
		StratifiedSampler* strat_sampler_jitter = new StratifiedSampler(std::ceil(pixel_indices / 1.0f), std::ceil(pixel_indices / 1.0f), true, seed);
		SobolSampler* sobol_sampler = new SobolSampler(pixel_indices, glm::ivec2(X_pixels * Y_pixels), RandomizeStrategy::FastOwen, seed);
		Sampler* current_sampler = indep_sampler;

		for (int x = 0; x < indep_grid.size(); x++)
		{
			for (int y = 0; y < indep_grid[x].size(); y++)
			{
				//at pixel(x,y) sample pixel_indices times
				glm::vec2 pixel(x, y);
				if (sampler_same_pixels)
					pixel = glm::vec2(0, 0);
				if(sampler_everypixel)
					current_sampler->StartPixelSample(pixel, 0, 0); //set index and dimension to 0 
				indep_grid[x][y].resize(pixel_indices);
				for (int s = 0; s < indep_grid[x][y].size(); s++)
				{
					if(!sampler_everypixel)
						current_sampler->StartPixelSample(pixel, (!sampler_same_pixels) ? s : 0, 0); //set index and dimension to 0 
					indep_grid[x][y][s] = current_sampler->Get2D();
				}
			}
		}

		//graph1 initialize
		glm::vec2 graph_scale(1850, 1850 / 1.5f);
		Graph2D graph1(graphics_manager, graph_scale.x, graph_scale.y);
		graph1.setbackgroundcolor(glm::vec3(1, 1, 1));
		graph1.setlinecolor(glm::vec3(0, 0, 255), 0);
		graph1.setlinecolor(glm::vec3(255, 0, 0), 1);
		graph1.setlinecolor(glm::vec3(0, 255, 0), 2);

		while (!glfwWindowShouldClose(graphics_manager.getWindow()))
		{
			//INPUT
			input_manager.HandleInput();
			glm::vec2 mouse_pos = graphics_manager.WindowToWorld(input_manager.getMousePosition());

			bool holding_not_first = (input_manager.getMouseButtonCurrent().x && !input_manager.getMouseButtonOnce().x);
			graph1.give_mouse_data(input_manager.getMousePosition(), input_manager.MouseScrollDirection(), holding_not_first, input_manager.getMouseButtonOnce().x);
	
			glm::vec3 graph_translate(0, 0, 0);
			graph1.setviewport(graph_translate.x, graph_translate.y, graph_scale.x, graph_scale.y);
			glm::vec3 graph_translate2(250, -250, 0);


			//graph predraw
			if (previous_method != current_method)
			{
				//if it uses graph
				if (current_method == 0 || current_method == 1 || current_method == 2 || current_method == 3 || current_method == 4
					|| current_method == 5 || current_method == 6 || current_method == 7 || current_method == 8)
				{
					graph1.clearalllines();
				}


				if (current_method == 0)
				{
					graph1.plotpoints(0, normal_samples);
					graph1.createfunction(1, -2, 2, .1, normal_dist);
				}
				else if (current_method == 1)
				{
					graph1.plotpoints(0, exp_sample_rejection);
					graph1.createfunction(1, 0, 2, .1, exponential_dist);
				}
				else if (current_method == 2)
				{
					for (int i = 0; i < geometric_dist.size(); i++)
					{
						float scale = 50.0f;
						if (i == 0)
							graph1.plotpoint(1, std::pair<float, float>(geom_xvalues[i], 0));
						else
							graph1.plotpoint(1, std::pair<float, float>(geom_xvalues[i], scale * geometric_dist[i - 1]));

						graph1.plotpoint(1, std::pair<float, float>(geom_xvalues[i], scale * geometric_dist[i]));

					}
					graph1.plotpoints(0, geometric_histogram);
				}
				else if (current_method == 3)
				{
					for (int i = 0; i < binomial_dist.size(); i++)
					{
						float scale = 50.0f;
						if (i == 0)
							graph1.plotpoint(1, std::pair<float, float>(binomial_xvalues[i], 0));
						else
							graph1.plotpoint(1, std::pair<float, float>(binomial_xvalues[i], scale * binomial_dist[i - 1]));

						graph1.plotpoint(1, std::pair<float, float>(binomial_xvalues[i], scale * binomial_dist[i]));

					}
					graph1.plotpoints(0, binomial_rejection_histogram);
				}
				if (current_method == 4)
				{
					graph1.plotpoints(0, linear_histogram);
					graph1.createfunction(1, 0, 1, .1, linear_func_normalized);
				}
				if (current_method == 5)
				{
					graph1.plotpoints(0, tent_histogram);
					graph1.createfunction(1,-tent_r, tent_r, .1, tent_func_normalized);
				}
				if (current_method == 6)
				{
					graph1.plotpoints(0, exp_histogram);
					graph1.createfunction(1, x_range_1, x_range_2, .1, exp_func_normalized);
				}
				if (current_method == 7)
				{
					graph1.plotpoints(0, gauss_histogram);
					graph1.createfunction(1, gauss_xrange_1, guass_xrange_2, .1, gauss_function_normalized);
				}
				if (current_method == 8)
				{
					graph1.plotpoints(0, logistic_histogram);
					graph1.createfunction(1, -5, 5, .1, logistics_function_normalized);
				}
			}

			graph1.DrawToTexture();

			//DRAW
			graphics_manager.Begin();

			//draw graph
			if (current_method == 0 || current_method == 1 || current_method == 2 || current_method == 3 || current_method == 4
				|| current_method == 5 || current_method == 6 || current_method == 7 || current_method == 8)
			{
				//draw the plotter texture
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), graph_translate);
				mat = glm::scale(mat, glm::vec3(graph_scale.x, graph_scale.y, 1)); //3:2 ratio
				DrawSubmit submit_1{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::TEXTURE, mat, glm::vec3(50,0,0) };
				submit_1.sampler_id = graph1.gettexture();
				graphics_manager.SubmitDraw(submit_1);

				glm::vec2 zzz = graph1.mousetographpos(input_manager.getMousePosition());
				std::string zzzz = std::to_string(zzz.x) + ", " + std::to_string(zzz.y);
				graphics_manager.RenderTextFont(zzzz, 450, 450, 0.5f, glm::vec3(0, 0, 0));
			}
			else if (current_method == 9)
			{
				float scaling_factor = 1000.0f;
				//draw a [0,1]^2 square
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0,0, 0));
				mat = glm::scale(mat, glm::vec3(scaling_factor, scaling_factor, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0),  });

				//plot all the sample points
				for (int i = 0; i < bilinear_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3((bilinear_samples[i].x-.5f)*scaling_factor, (bilinear_samples[i].y - .5f) *scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 10)
			{
				float scaling_factor = 500.0f;
				//draw disk
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(scaling_factor, scaling_factor, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0) });
				
				//plot all the sample points
				for (int i = 0; i < diskpolar_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(diskpolar_samples[i].x * scaling_factor, diskpolar_samples[i].y * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 11)
			{
				float scaling_factor = 500.0f;
				//draw disk
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(scaling_factor, scaling_factor, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0) });

				//plot all the sample points
				for (int i = 0; i < diskconcentric_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(diskconcentric_samples[i].x * scaling_factor, diskconcentric_samples[i].y * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 12)
			{
				float scaling_factor = 500.0f;
				//draw disk
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(scaling_factor, scaling_factor, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0) });

				//plot all the sample points
				for (int i = 0; i < hemisphere_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(hemisphere_samples[i].x * scaling_factor, hemisphere_samples[i].y * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 13)
			{
				float scaling_factor = 250.0f;
				//draw disk
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(250, 0, 0));
				mat = glm::scale(mat, glm::vec3(scaling_factor, scaling_factor, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0) });

				mat = glm::translate(glm::mat4(1.0f), glm::vec3(-250, 0, 0));
				mat = glm::scale(mat, glm::vec3(scaling_factor, scaling_factor, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0) });

				//plot all the sample points
				for (int i = 0; i < sphere_samples.size(); i++)
				{
					glm::mat4 mat(1.0f);
					if (sphere_samples[i].z <= 0)
						mat = glm::translate(glm::mat4(1.0f), glm::vec3(sphere_samples[i].x* scaling_factor - 250.0f, sphere_samples[i].y * scaling_factor, 0));
					else
						mat = glm::translate(glm::mat4(1.0f), glm::vec3(sphere_samples[i].x * scaling_factor + 250, sphere_samples[i].y * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 14)
			{
				float scaling_factor = 500.0f;
				//draw disk
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(scaling_factor, scaling_factor, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0) });

				//plot all the sample points
				for (int i = 0; i < cosinhemisphere_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(cosinhemisphere_samples[i].x * scaling_factor, cosinhemisphere_samples[i].y * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 15)
			{
				float scaling_factor = 500.0f;
				//draw disk
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(scaling_factor, scaling_factor, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0) });

				//plot all the sample points
				for (int i = 0; i < cone_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(cone_samples[i].x * scaling_factor, cone_samples[i].y * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 16)
			{
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(sampler_scaling, sampler_scaling, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0) });

				float starting_x = -.5 * sampler_scaling;
				float ending_x = .5 * sampler_scaling;
				float starting_y = .5 * sampler_scaling;
				float ending_y = -.5 * sampler_scaling;
				float delta_x = (1.0f / (float)X_pixels) * sampler_scaling;
				float delta_y = (1.0f / (float)Y_pixels) * sampler_scaling;
				float uniform_delta_x = 1.0f / (float)X_pixels;
				float uniform_delta_y = 1.0f / (float)Y_pixels;
				for (int x = 0; x < indep_grid.size(); x++)
				{
					float current_x = starting_x + delta_x * x;
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(255,255,255), 
						glm::vec4(current_x,starting_y,current_x,ending_y) });
					for (int y = 0; y < indep_grid[x].size(); y++)
					{
						float current_y = starting_y - delta_y * y;
						graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::LINE, SHADER_TYPE::MAIN, glm::mat4(1.0f), glm::vec3(255,255,255),
						glm::vec4(starting_x,current_y,ending_x,current_y) });
						for (int s = 0; s < indep_grid[x][y].size(); s++)
						{
							//map [0,1]^2 to the [0,1]^2 grid spot
							//map to the correct grid location, shrink to the size of 1 grid, then move it to bottomleft of specific grid
							glm::vec2 uniform_pos = indep_grid[x][y][s];
							glm::vec2 shrink_pos = uniform_pos * glm::vec2(uniform_delta_x, uniform_delta_y);
							glm::vec2 pos = shrink_pos + glm::vec2(0 + uniform_delta_x * x, 1 - uniform_delta_y * (y + 1));

							//color coat based off s
							glm::vec3 color = glm::vec3(255, 0, 0);
							if(sampler_colorcoat)
								color = mapindextocolor(s);

							glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3((pos.x-.5f) * sampler_scaling, (pos.y-.5f) * sampler_scaling, 0));
							mat = glm::scale(mat, glm::vec3(sampler_point_scaling, sampler_point_scaling, 1));
							graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, color });
						}
					}
				}
			}

			graphics_manager.Render();


			//IMGUI
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();

			ImGui::NewFrame();
			//ImGui::ShowDemoWindow();
			ImGui::Begin("Sampler Test Solver app");
			
			
			previous_method = current_method;
			ImGui::Combo(name, &current_method, sampling_methods.data(), sampling_methods.size());
			
			ImGui::Text("---Uniform Pixel Sampler config---");
			bool generate_new_samples = false;
			ImGui::DragFloat("square scaling", &sampler_scaling, 1.0f, 1.0f, 1500.0f);
			ImGui::DragFloat("point scaling", &sampler_point_scaling, 1.0f, 1.0f, 15.0f);
			ImGui::Checkbox("color coat", &sampler_colorcoat);
			if (ImGui::Checkbox("reset every pixel instead of pixel index", &sampler_everypixel))
				generate_new_samples = true;
			if (ImGui::Checkbox("same pixel samples", &sampler_same_pixels))
				generate_new_samples = true;

			if (ImGui::Combo("sampler_type", &sampler_type, sampler_types.data(), sampler_types.size()))
			{
				switch (sampler_type)
				{
				case 0: current_sampler = indep_sampler; break;
				case 1: current_sampler = strat_sampler; break;
				case 2: current_sampler = strat_sampler_jitter; break;
				case 3: current_sampler = sobol_sampler; break;
				default: break;
				}
				generate_new_samples = true;
			}

			if (ImGui::Button("randomize"))
			{
				seed = Helper::GetRandomNumber(0, 1000000.0);
				generate_new_samples = true;
				
				delete indep_sampler;
				delete strat_sampler;
				delete strat_sampler_jitter;
				delete sobol_sampler;
				indep_sampler = new IndependentSampler(pixel_indices, seed);
				strat_sampler = new StratifiedSampler(std::ceil(pixel_indices / 1.0f), std::ceil(pixel_indices / 1.0f), false, seed);
				strat_sampler_jitter = new StratifiedSampler(std::ceil(pixel_indices / 1.0f), std::ceil(pixel_indices / 1.0f), true, seed);
				sobol_sampler = new SobolSampler(pixel_indices, glm::ivec2(X_pixels * Y_pixels), RandomizeStrategy::PermuteDigits, seed);
				switch (sampler_type)
				{
				case 0: current_sampler = indep_sampler; break;
				case 1: current_sampler = strat_sampler; break;
				case 2: current_sampler = strat_sampler_jitter; break;
				case 3: current_sampler = sobol_sampler; break;
				default: break;
				}
			}

			int previous_pixel_index = pixel_indices;
			if (ImGui::InputInt("per pixel samples", &pixel_indices, 1)) 
			{
				generate_new_samples = true;
				//stratified gets stuck on has if we go over
				if (sampler_type == 1 && pixel_indices >= current_sampler->SamplesPerPixel())
				{
					pixel_indices = previous_pixel_index;
				}
			}

			if (generate_new_samples)
			{
				for (int x = 0; x < indep_grid.size(); x++)
				{
					for (int y = 0; y < indep_grid[x].size(); y++)
					{
						//at pixel(x,y) sample pixel_indices times
						glm::vec2 pixel(x, y);
						if (sampler_same_pixels)
							pixel = glm::vec2(0, 0);
						if (sampler_everypixel)
							current_sampler->StartPixelSample(pixel, 0, 0); //set index and dimension to 0 
						indep_grid[x][y].clear();
						indep_grid[x][y].resize(pixel_indices);
						for (int s = 0; s < indep_grid[x][y].size(); s++)
						{
							if (!sampler_everypixel)
								current_sampler->StartPixelSample(pixel, (!sampler_same_pixels) ? s : 0, 0); //set index and dimension to 0 
							indep_grid[x][y][s] = current_sampler->Get2D();
							//std::cout << indep_grid[x][y][s].x << " " << indep_grid[x][y][s].y << std::endl;
						}
					}
				}
			}

			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);


			ImGui::End();

			ImGui::Render();
			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

			graphics_manager.End();
		}

		delete indep_sampler;
		delete strat_sampler;
		delete strat_sampler_jitter;
	}

	glm::vec3 mapindextocolor(int index)
	{
		glm::vec3 color(0, 0, 0);
		index = index % 13;
		switch (index)
		{
		case 0: color = glm::vec3(255, 0, 0); break;
		case 1: color = glm::vec3(0, 255, 255); break;
		case 2: color = glm::vec3(255, 255, 0); break;
		case 3: color = glm::vec3(255, 0, 255); break;
		case 4: color = glm::vec3(0, 255, 0); break;
		case 5: color = glm::vec3(0, 0, 255); break;
		case 6: color = glm::vec3(0, 51, 51); break;
		case 7: color = glm::vec3(200, 200, 255); break;
		case 8: color = glm::vec3(100, 100, 100); break;
		case 9: color = glm::vec3(102, 0, 0); break;
		case 10: color = glm::vec3(200, 255, 150); break;
		case 11: color = glm::vec3(0, 153, 0); break;
		case 12: color = glm::vec3(0, 0, 153); break;
		default: break;
		}

		return color;
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