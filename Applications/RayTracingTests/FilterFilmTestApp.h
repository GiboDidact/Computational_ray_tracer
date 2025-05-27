#pragma once
#include "../../Graphics/GraphicsManager.h"
#include "../../Graphics/Graph2D.h"
#include "../../RayTracer/Sampling.h"
#include "../../ThirdParty/pbrv4/samplers.h"
#include "../../ThirdParty/pbrv4/filters.h"
#include "../../ThirdParty/pbrv4/pixelsensor.h"
#include "../../ThirdParty/pbrv4/color.h"

using namespace pbrt;

class FilterFilmTestApp
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		FilterFilmTestApp* handler = reinterpret_cast<FilterFilmTestApp*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<FilterFilmTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<FilterFilmTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<FilterFilmTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<FilterFilmTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons;
		if (graphics_manager.InitWindow(1850, 1850 / 1.5f, "FilterFilmTest App", icons) == EXIT_FAILURE)
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
		std::vector<const char*> sampling_methods = { "box_filter_x", "box_filter_y", "box_filter_samples",
		"tri_filter_x", "tri_filter_samples", "gauss_filter_x", "guass_filter_samples", "lanczos_filter_x", "lanczos_filter_samples",
		"visible wavelength graph", "visible wavelength samples", "rgbsigmoid", "sigmoid_s", "sigmoid_polynomial", "sRGB_to_spectral"
		, "albedo_spectral" };
		int current_method = 14;
		int previous_method = -1;
		const char* name = "sampling methods";

		//FILTERS
		//for each filter: graph it, sample and draw histogram, then display radius and integrand
		//These are 2D but seperable, so do X and Y seperately

		//BOX
		float box_r = 5.0f;
		BoxFilter box_filter(glm::vec2(box_r, box_r));
		auto box_filter_x = [=](float x) -> float { return box_filter.Evaluate(glm::vec2(x, 0)); };
		auto box_filter_y = [=](float y) -> float { return box_filter.Evaluate(glm::vec2(0, y)); };
		std::cout << box_filter.ToString() << " integral: " << box_filter.Integral() << std::endl;
		std::vector<glm::vec2> box_filter_samples;
		box_filter_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			FilterSample fs = box_filter.Sample(glm::vec2(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1)));
			box_filter_samples.push_back(fs.p);
		}

		//TRI
		float tri_r = 3.0f;
		TriangleFilter tri_filter(glm::vec2(tri_r,tri_r));
		auto tri_filter_x = [=](float x) -> float { return tri_filter.Evaluate(glm::vec2(x, 0)); };
		auto tri_filter_y = [=](float y) -> float { return tri_filter.Evaluate(glm::vec2(0, y)); };
		std::cout << tri_filter.ToString() << " integral: " << tri_filter.Integral() << std::endl;
		std::vector<glm::vec2> tri_filter_samples;
		tri_filter_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			FilterSample fs = tri_filter.Sample(glm::vec2(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1)));
			tri_filter_samples.push_back(fs.p);
		}

		//GAUSS
		float gauss_r = 3.0f;
		GaussianFilter gauss_filter(glm::vec2(tri_r, tri_r), 0.5f);
		auto gauss_filter_x = [=](float x) -> float { return gauss_filter.Evaluate(glm::vec2(x, 0)); };
		auto gauss_filter_y = [=](float y) -> float { return gauss_filter.Evaluate(glm::vec2(0, y)); };
		std::cout << gauss_filter.ToString() << " integral: " << gauss_filter.Integral() << std::endl;
		std::vector<glm::vec2> gauss_filter_samples;
		gauss_filter_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			FilterSample fs = gauss_filter.Sample(glm::vec2(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1)));
			gauss_filter_samples.push_back(fs.p);
		}

		//LANCZOS
		float lancz_r = 3.0f;
		LanczosSincFilter lancz_filter(glm::vec2(tri_r, tri_r));
		auto lancz_filter_x = [=](float x) -> float { return lancz_filter.Evaluate(glm::vec2(x, 0)); };
		auto lancz_filter_y = [=](float y) -> float { return lancz_filter.Evaluate(glm::vec2(0, y)); };
		std::cout << lancz_filter.ToString() << " integral: " << lancz_filter.Integral() << std::endl;
		std::vector<glm::vec2> lancz_filter_samples;
		lancz_filter_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			FilterSample fs = lancz_filter.Sample(glm::vec2(Helper::GetRandomNumber(0, 1), Helper::GetRandomNumber(0, 1)));
			lancz_filter_samples.push_back(fs.p);
		}

		//PIXEL SENSOR
		//exposure, spectral/color conversion, whitebalancing
		Spectra::Init();
		RGBToSpectrumTable::Init();
		RGBColorSpace::Init();

		RGBColorSpace _srgb(glm::vec2(.64, .33), glm::vec2(.3, .6), glm::vec2(.15, .06),
			GetNamedSpectrum("stdillum-D65"), nullptr);
		Spectrum* illumA = GetNamedSpectrum("stdillum-A");
		Spectrum* illumD = GetNamedSpectrum("stdillum-D65");
		Spectrum* canon_r = GetNamedSpectrum("canon_eos_100d_r");
		Spectrum* canon_g = GetNamedSpectrum("canon_eos_100d_g");
		Spectrum* canon_b = GetNamedSpectrum("canon_eos_100d_b");
		Spectrum* sony_r = GetNamedSpectrum("sony_ilce_9_r");
		Spectrum* sony_g = GetNamedSpectrum("sony_ilce_9_g");
		Spectrum* sony_b = GetNamedSpectrum("sony_ilce_9_b");

		SampledWavelengths wave_lengths_2 = SampledWavelengths::SampleVisible(Helper::GetRandomNumber(0, 1));
		SampledSpectrum spectrum_2 = illumA->Sample(wave_lengths_2);// wave_lengths_2.ToSampled(illumA);
		XYZ spec_2_xyz = spectrum_2.ToXYZ(wave_lengths_2);
		//RGB spec_2_rgb = spectrum_2.ToRGB(wave_lengths_2);
		
		PixelSensor sensor_xyz(&_srgb, illumD, 1.0f/ CIE_Y_integral);
		RGB sensor_rgb = sensor_xyz.ToSensorRGB(spectrum_2, wave_lengths_2);
		//std::cout << spec_2_xyz.ToString() << " " << sensor_rgb.ToString() << std::endl;

		PixelSensor sensor_canon(canon_r, canon_g, canon_b, &_srgb, illumA, 1.0f / CIE_Y_integral);
		RGB canon_rgb = sensor_canon.ToSensorRGB(spectrum_2, wave_lengths_2);
		XYZ canon_xyz = sensor_canon.XYZFromSensorRGB * canon_rgb.getglm();
		//std::cout << canon_rgb.ToString() << std::endl;
		std::cout << canon_xyz.ToString() << std::endl;

		PixelSensor sensor_sony(sony_r, sony_g, sony_b, &_srgb, illumA, 1.0f / CIE_Y_integral);
		RGB sony_rgb = sensor_sony.ToSensorRGB(spectrum_2, wave_lengths_2);
		XYZ sony_xyz = sensor_sony.XYZFromSensorRGB * sony_rgb.getglm();
		std::cout << sony_xyz.ToString() << std::endl;

		//show the sensor response curve sampling functions
		auto visible_pdf = [=](float x) -> float { return 10000*VisibleWavelengthsPDF(x+360); };
		std::vector<float> visible_samples;
		visible_samples.reserve(1000);
		for (int i = 0; i < 10000; i++)
		{
			float val = SampleVisibleWavelengths(Helper::GetRandomNumber(0, 1));
			visible_samples.push_back(val);
		}
		std::vector<float> visible_histogram = generateContinuousSampleHistogram(visible_samples, 360, 830, 100, 1500.0f, -1);

		//we have the estimator for filtering integral


		//RGB TO SPECTRUM
		float c0 = .0001, c1 = .00005, c2 = .0001;
		RGBSigmoidPolynomial sigmoid(c0,c1,c2);
		std::function<float(float)> sigmoid_func = [=](float x) -> float { return 1 * sigmoid(x + 360); };
		std::function<float(float)> s_func = [=](float x) -> float { return 1 * sigmoid.s(x); };
		std::function<float(float)> sigpoly_func = [=](float x) -> float { return .1 * sigmoid.polynomial(x); };

		//match (.7,.5,.8) purple, (.25,.44,.33) green, (.36,.275,.21) brown
		//try uniform rgb spectral should be constant
		//map rgb to spectral and spectral back, they should be the same
		RGBSigmoidPolynomial red_poly = RGBToSpectrumTable::sRGB->operator()(RGB(.7, .5, .8));
		std::function<float(float)> tospectral_func = [=](float x) -> float { return 100 * red_poly(x); };


		RGBAlbedoSpectrum* rgb_to_sprectrum = new RGBAlbedoSpectrum(*RGBColorSpace::sRGB, RGB(1, 0, 0));
		std::function<float(float)> albedo_func = [=](float x) -> float { return 1* rgb_to_sprectrum->Query(x); }; 
		
		SampledSpectrum albed_samples = rgb_to_sprectrum->Sample(wave_lengths_2);// wave_lengths_2.ToSampled(rgb_to_sprectrum);
		RGB val = albed_samples.ToRGB(wave_lengths_2, *RGBColorSpace::sRGB);
		std::cout << val.ToString() << std::endl;

		//graph1 initialize
		glm::vec2 graph_scale(1850, 1850 / 1.5f);
		Graph2D graph1(graphics_manager, graph_scale.x, graph_scale.y);
		graph1.setbackgroundcolor(glm::vec3(1, 1, 1));
		graph1.setlinecolor(glm::vec3(0, 0, 255), 0);
		graph1.setlinecolor(glm::vec3(255, 0, 0), 1);
		graph1.setlinecolor(glm::vec3(0, 255, 0), 2);


		float tri_filter_yvalue = 0.0f;
		float gauss_filter_yvalue = 0.0f;
		float lanczos_filter_yvalue = 0.0f;
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
					|| current_method == 5 || current_method == 6 || current_method == 7 || current_method == 8 || current_method == 9 || current_method == 11
					|| current_method == 12 || current_method == 13 || current_method == 14 || current_method == 15)
				{
					graph1.clearalllines();
				}

				if (current_method == 0)
				{
					//graph1.plotpoints(0, normal_samples);
					graph1.createfunction(1, -box_r, box_r, 1, box_filter_x);
				}
				else if (current_method == 1)
				{
					//graph1.plotpoints(0, exp_sample_rejection);
					graph1.createfunction(1, -box_r, box_r, 1, box_filter_y);
				}
				else if (current_method == 3)
				{
					graph1.createfunction(1, -tri_r, tri_r, 1, tri_filter_x);
				}
				else if (current_method == 5)
				{
					graph1.createfunction(1, -gauss_r, gauss_r, .01, gauss_filter_x);
				}
				else if (current_method == 7)
				{
					graph1.createfunction(1, -lancz_r, lancz_r, .01, lancz_filter_x);
				}
				else if (current_method == 9)
				{
					graph1.plotpoints(0, visible_histogram);
					graph1.createfunction(1, 0, 840-360, 1, visible_pdf);
				}
				else if (current_method == 11)
				{
					graph1.createfunction(1, 0, 840 - 360, 1, sigmoid_func);
				}
				else if (current_method == 12)
				{
					graph1.createfunction(1, -10, 10, 1, s_func);
				}
				else if (current_method == 13)
				{
					graph1.createfunction(1, 360, 830, 1, sigpoly_func);
				}
				else if (current_method == 14)
				{
					graph1.createfunction(1, 360, 830, 1, tospectral_func);
				}
				else if (current_method == 15)
				{
					graph1.createfunction(1, 360, 830, 1, albedo_func);
				}
			}

			graph1.DrawToTexture();

			//DRAW
			graphics_manager.Begin();

			//draw graph
			if (current_method == 0 || current_method == 1 || current_method == 3 || current_method == 5 || current_method == 7 || current_method == 9\
				|| current_method == 11 || current_method == 12 || current_method == 13 || current_method == 14 || current_method == 15)
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
			else if (current_method == 2)
			{
				float scaling_factor = 1000.0f;
				//draw a [0,1]^2 square
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(1000.0f, 1000.0f, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0), });

				//plot all the sample points
				for (int i = 0; i < box_filter_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3((box_filter_samples[i].x / (2 * box_r)) * scaling_factor, (box_filter_samples[i].y / (2 * box_r)) * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 4)
			{
				float scaling_factor = 1000.0f;
				//draw a [0,1]^2 square
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(1000.0f, 1000.0f, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0), });

				//plot all the sample points
				for (int i = 0; i < tri_filter_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3((tri_filter_samples[i].x / (2 * tri_r)) * scaling_factor, (tri_filter_samples[i].y / (2 * tri_r)) * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 6)
			{
				float scaling_factor = 1000.0f;
				//draw a [0,1]^2 square
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(1000.0f, 1000.0f, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0), });

				//plot all the sample points
				for (int i = 0; i < gauss_filter_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3((gauss_filter_samples[i].x / (2 * gauss_r)) * scaling_factor, (gauss_filter_samples[i].y / (2 * gauss_r)) * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 8)
			{
				float scaling_factor = 1000.0f;
				//draw a [0,1]^2 square
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(1000.0f, 1000.0f, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0), });

				//plot all the sample points
				for (int i = 0; i < lancz_filter_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3((lancz_filter_samples[i].x / (2 * lancz_r)) * scaling_factor, (lancz_filter_samples[i].y / (2 * lancz_r)) * scaling_factor, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}
			else if (current_method == 10)
			{
				float scaling_factor = 1000.0f;
				//draw a [0,1]^2 square
				glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
				mat = glm::scale(mat, glm::vec3(1000.0f, 1000.0f, 1));
				graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::MAIN, mat, glm::vec3(0,0,0), });

				//plot all the sample points
				for (int i = 0; i < visible_samples.size(); i++)
				{
					glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(((visible_samples[i]-595.0) / (2* -235.0f)) * scaling_factor, 0, 0));
					mat = glm::scale(mat, glm::vec3(1, 1, 1));
					graphics_manager.SubmitDraw(DrawSubmit{ GEOMETRY_TYPE::CIRCLE, SHADER_TYPE::MAIN, mat, glm::vec3(255,0,0) });
				}
			}

			graphics_manager.Render();


			//IMGUI
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();

			ImGui::NewFrame();
			ImGui::ShowDemoWindow();
			ImGui::Begin("FilterFilmTest app");


			previous_method = current_method;
			ImGui::Combo(name, &current_method, sampling_methods.data(), sampling_methods.size());

			if (ImGui::SliderFloat("tri_filter_y_value", &tri_filter_yvalue, -tri_r, tri_r))
			{
				if (current_method == 3)
				{
					auto tri_filter_x = [=](float x) -> float { return tri_filter.Evaluate(glm::vec2(x, tri_filter_yvalue)); };
					graph1.createfunction(1, -tri_r, tri_r, 1, tri_filter_x);
				}
			}
			if (ImGui::SliderFloat("gauss_filter_y_value", &gauss_filter_yvalue, -gauss_r, gauss_r))
			{
				if (current_method == 5)
				{
					auto gaus_filter_x2 = [=](float x) -> float { return gauss_filter.Evaluate(glm::vec2(x, gauss_filter_yvalue)); };
					graph1.createfunction(1, -gauss_r, gauss_r, .1, gaus_filter_x2);
				}
			}
			if (ImGui::SliderFloat("lanczos_filter_y_value", &lanczos_filter_yvalue, -lancz_r, lancz_r))
			{
				if (current_method == 7)
				{
					auto lanczos_filter_x2 = [=](float x) -> float { return lancz_filter.Evaluate(glm::vec2(x, lanczos_filter_yvalue)); };
					graph1.createfunction(1, -lancz_r, lancz_r, .1, lanczos_filter_x2);
				}
			}

			if (ImGui::SliderFloat("c0", &c0, -10, 10) || ImGui::SliderFloat("c1", &c1, -10, 10) || ImGui::SliderFloat("c2", &c2, -10, 10))
			{
				if (current_method == 11)
				{
					RGBSigmoidPolynomial sigmoid2(c0, c1, c2);
					std::cout << sigmoid2(360) << " " << sigmoid2(500) << " " << sigmoid2(800) << std::endl;
					sigmoid_func = [=](float x) -> float { return 1000 * sigmoid2(x + 360); };
				};
			}

			float colors[3];
			if (ImGui::ColorPicker3("rgb to spectral", colors))
			{
				if (current_method == 14)
				{
					red_poly = RGBToSpectrumTable::sRGB->operator()(RGB(colors[0], colors[1], colors[2]));
					std::cout << red_poly.MaxValue() <<" "<<red_poly(500)<< std::endl;
					tospectral_func = [=](float x) -> float { return 100 * red_poly(x); };
					graph1.createfunction(1, 360, 830, 1, tospectral_func);
				}
				else if (current_method == 15)
				{
					delete rgb_to_sprectrum;

					rgb_to_sprectrum = new RGBAlbedoSpectrum(*RGBColorSpace::sRGB, RGB(colors[0], colors[1], colors[2]));
					RGBUnboundedSpectrum* unbounded = new RGBUnboundedSpectrum(*RGBColorSpace::sRGB, RGB(colors[0], colors[1], colors[2]));
					std::cout << rgb_to_sprectrum->MaxValue()<< std::endl;
					albedo_func = [=](float x) -> float { return 100 * rgb_to_sprectrum->Query(x); };
					graph1.createfunction(1, 360, 830, 1, albedo_func);

					XYZ xyz_1 = RGBColorSpace::sRGB->ToXYZ(RGB(colors[0], colors[1], colors[2]));
					XYZ val_xyz = SpectrumToXYZ(rgb_to_sprectrum);
					RGB val = RGBColorSpace::sRGB->ToRGB(val_xyz);
					SampledSpectrum albed_samples = rgb_to_sprectrum->Sample(wave_lengths_2);// wave_lengths_2.ToSampled(rgb_to_sprectrum);
					RGB val2 = albed_samples.ToRGB(wave_lengths_2, *RGBColorSpace::sRGB);
					std::cout <<"full spectral-->srgb "<< val.ToString() << std::endl;
					std::cout <<"sampled spectral-->srgb "<< val2.ToString() << std::endl;
					std::cout << RGB(colors[0], colors[1], colors[2]).ToString() << std::endl;
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