#pragma once
#include "../../Graphics/GraphicsManager.h"
#include "../../ThirdParty/pbrv4/spectrum.h"
#include "../../ThirdParty/pbrv4/colorspace.h"
#include "../../ThirdParty/pbrv4/samplers.h"
#include "../../Graphics/Graph2D.h"
using namespace pbrt;

class SpectrumColorTestApp
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		SpectrumColorTestApp* handler = reinterpret_cast<SpectrumColorTestApp*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<SpectrumColorTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<SpectrumColorTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<SpectrumColorTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<SpectrumColorTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons;
		if (graphics_manager.InitWindow(1850, 1850 / 1.5f, "SpectrumColorTest App", icons) == EXIT_FAILURE)
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
		std::cout << "\n------SPECTRUM GRAPHS-------\n";
		Spectra::Init();
	
		glm::vec2 graph_scale(1850, 1850 / 1.5f);
		Graph2D graph1(graphics_manager, graph_scale.x, graph_scale.y);
		graph1.setbackgroundcolor(glm::vec3(1, 1, 1));
		graph1.setlinecolor(glm::vec3(0, 0, 255), 0);
		graph1.setlinecolor(glm::vec3(255, 0, 0), 1);
		graph1.setlinecolor(glm::vec3(0, 255, 0), 2);
		graph1.setlinecolor(glm::vec3(0, 255, 255), 3);
		graph1.setlinecolor(glm::vec3(0, 0, 0), 4);
		graph1.move_to_pos(glm::vec2(550, 0));
		graph1.scale_cam(glm::vec2(Lambda_max - Lambda_min, (Lambda_max - Lambda_min)/1.5f));

		//Spectrums
		float T = 3000.0f;
		BlackbodySpectrum body_spec(T);
		std::cout << "blackbody max: " << body_spec.MaxValue() << std::endl;
		auto plank_func = [=](float lambda) -> float { return 500*body_spec.Query(lambda); };
		//graph1.createfunction(0, Lambda_min, Lambda_max, .1, plank_func);

		ConstantSpectrum const_spect(1);
		std::cout << "cont max: " << const_spect.MaxValue() << std::endl;
		auto const_func = [=](float lambda) -> float { return const_spect.Query(lambda); };
		//graph1.createfunction(1, Lambda_min, Lambda_max, .1, const_func);

		auto linear_func = [=](float x) -> float { return .1*(x-Lambda_min); };
		DenselySampledSpectrum dense_spect = DenselySampledSpectrum::SampleFunction(linear_func);
		auto dense_func = [=](float lambda) -> float { return dense_spect.Query(lambda); };
		//graph1.createfunction(2, Lambda_min, Lambda_max, .1, dense_func);

		DenselySampledSpectrum dense_spect2(&dense_spect);
		dense_spect2.Scale(2);
		auto dense_func2 = [=](float lambda) -> float { return dense_spect2.Query(lambda); };
		//graph1.createfunction(3, Lambda_min, Lambda_max, .1, dense_func2);

		std::vector<float> lambda_1 = { 360.0, 400.0, 600.0, 800.0, 860.0 };
		std::vector<float> values_1 = { 10.0, 50.0, 10.0, 40.0, 80.0 };
		std::span<const float> lol(lambda_1);
		std::span<const float> lol2(values_1);
		PiecewiseLinearSpectrum piece_spect(lol, lol2);
		piece_spect.Scale(2);
		std::cout << "piece max: " << piece_spect.MaxValue() << std::endl;
		auto piece_func = [=](float lambda) -> float { return piece_spect.Query(lambda); };
		//graph1.createfunction(4, Lambda_min, Lambda_max, .1, piece_func);

		//spectrum data tables: X,Y,Z, light emission, glass, camera sensitivity
		const Spectrum* X = &Spectra::X();
		const Spectrum* Y = &Spectra::Y();
		const Spectrum* Z = &Spectra::Z();

		auto X_func = [=](float lambda) -> float { return 500 * X->Query(lambda); };
		auto Y_func = [=](float lambda) -> float { return 500 * Y->Query(lambda); };
		auto Z_func = [=](float lambda) -> float { return 500 * Z->Query(lambda); };

		//graph1.clearalllines();
		//graph1.createfunction(0, Lambda_min, Lambda_max, .1, X_func);
		//graph1.createfunction(1, Lambda_min, Lambda_max, .1, Y_func);
		//graph1.createfunction(2, Lambda_min, Lambda_max, .1, Z_func);

		Spectrum* illumA = GetNamedSpectrum("stdillum-A");
		Spectrum* illumD = GetNamedSpectrum("stdillum-D65");
		Spectrum* illumF = GetNamedSpectrum("stdillum-F9");
		Spectrum* glassBK7 = GetNamedSpectrum("glass-BK7");
		Spectrum* metalAG  = GetNamedSpectrum("metal-Ag-eta");
		Spectrum* canon_r  = GetNamedSpectrum("canon_eos_100d_r");
		Spectrum* canon_g  = GetNamedSpectrum("canon_eos_100d_g");
		Spectrum* canon_b  = GetNamedSpectrum("canon_eos_100d_b"); 

		Spectrum* sony_r = GetNamedSpectrum("sony_ilce_9_r");
		Spectrum* sony_g = GetNamedSpectrum("sony_ilce_9_g");
		Spectrum* sony_b = GetNamedSpectrum("sony_ilce_9_b");
	
		std::vector<float> grasslambda = { 360.0, 500.0, 550.0, 690.0, 750.0, 830 };
		std::vector<float> grassvalues = { 3.5, 5, 11.0, 5.6, 46.0, 50.0};
		std::span<const float> glol(grasslambda);
		std::span<const float> glol2(grassvalues);
		PiecewiseLinearSpectrum grass_spec(glol, glol2);
		auto grass_func = [=](float lambda) -> float { return 50 * grass_spec.Query(lambda); };


		auto illumA_func = [=](float lambda) -> float { return 500 * illumA->Query(lambda); };
		auto illumD_func = [=](float lambda) -> float { return 500 * illumD->Query(lambda); };
		auto illumF_func = [=](float lambda) -> float { return 500 * illumF->Query(lambda); };
		auto glassBK7_func = [=](float lambda) -> float { return 500 * glassBK7->Query(lambda); };
		auto metalAG_func = [=](float lambda) -> float { return 500 * metalAG->Query(lambda); };
		auto canonr_func = [=](float lambda) -> float { return 500 * canon_r->Query(lambda); };
		auto canong_func = [=](float lambda) -> float { return 500 * canon_g->Query(lambda); };
		auto canonb_func = [=](float lambda) -> float { return 500 * canon_b->Query(lambda); };

		auto sonyr_func = [=](float lambda) -> float { return 500 * sony_r->Query(lambda); };
		auto sonyg_func = [=](float lambda) -> float { return 500 * sony_g->Query(lambda); };
		auto sonyb_func = [=](float lambda) -> float { return 500 * sony_b->Query(lambda); };


		graph1.clearalllines();
		graph1.createfunction(0, Lambda_min, Lambda_max, .1, plank_func);
		graph1.createfunction(1, Lambda_min, Lambda_max, .1, const_func);
		graph1.createfunction(2, Lambda_min, Lambda_max, .1, dense_func);
		graph1.createfunction(3, Lambda_min, Lambda_max, .1, dense_func2);
		graph1.createfunction(4, Lambda_min, Lambda_max, .1, piece_func);

		graph1.clearalllines();
		//graph1.createfunction(0, Lambda_min, Lambda_max, .1, illumD_func);
		//graph1.createfunction(1, Lambda_min, Lambda_max, .1, illumA_func);
		//graph1.createfunction(2, Lambda_min, Lambda_max, .1, illumF_func);

		//graph1.clearalllines();
		//graph1.createfunction(0, Lambda_min, Lambda_max, .1, X_func);
		//graph1.createfunction(1, Lambda_min, Lambda_max, .1, Y_func);
		//graph1.createfunction(2, Lambda_min, Lambda_max, .1, Z_func);
		

		//SAMPLING WAVELENGTH
		std::cout << "\n----------SAMPLING WAVELENGTH----------\n";
		StratifiedSampler strat_sampler(1850, 1850 / 1.5f, false, 1);
		strat_sampler.StartPixelSample(glm::vec2(0, 0), 0, 0);
		SampledWavelengths wave_lengths_1 = SampledWavelengths::SampleUniform(strat_sampler.Get1D());
		std::cout << wave_lengths_1.ToString() << std::endl;

		SampledWavelengths wave_lengths_2 = SampledWavelengths::SampleVisible(Helper::GetRandomNumber(0,1));
		//wave_lengths_2.TerminateSecondary();
		std::cout << wave_lengths_2.ToString() << std::endl;

		//SPECTRUM CONVERSION XYZ/PHOTOMETRIC
		std::cout << "------SPECTRUM CONVERSION------\n";
		SampledSpectrum spectrum_2 = illumA->Sample(wave_lengths_2);// wave_lengths_2.ToSampled(illumA);
		XYZ xyz = spectrum_2.ToXYZ(wave_lengths_2);
		//RGBColorSpace space_1;
		//RGB rgb = spectrum_2.ToRGB(wave_lengths_2, space_1);
		float y = spectrum_2.y(wave_lengths_2);
		std::cout << "xyz: " << xyz.ToString() << " y: " << y << std::endl;
		XYZ xyz_full = SpectrumToXYZ(illumA);
		std::cout << "xyz: " << xyz_full.ToString() << std::endl;
		float photo_1 = SpectrumToPhotometric(illumA);
		float photo_2 = SpectrumToPhotometric(illumF);
		std::cout << "lumA: " << photo_1 << " lumD: " << photo_2 << std::endl;


		//COLOR SPACE
		std::cout << "-------COLOR SPACE------\n";
		RGBColorSpace::Init();
		RGBColorSpace::sRGB;
		RGBColorSpace::GetNamed("srgb");
		//RGBColorSpace::Lookup();
		//ConvertRGBColorSpace();
	
		RGBColorSpace _srgb(glm::vec2(.64, .33), glm::vec2(.3, .6), glm::vec2(.15, .06),
			GetNamedSpectrum("stdillum-D65"), nullptr);
		RGBColorSpace _DCI_P3(RGBColorSpace(glm::vec2(.68, .32), glm::vec2(.265, .690), glm::vec2(.15, .06),
			GetNamedSpectrum("stdillum-D65"), nullptr));

		XYZ xyz_srgb = _srgb.ToXYZ(RGB(.2, .5, 1));
		RGB srgb_xyz = _srgb.ToRGB(XYZ(.5, .4, .5));
		std::cout << "srgb to xyz: " << xyz_srgb.ToString() << std::endl;
		std::cout << "xyz to srgb: " << srgb_xyz.ToString() << std::endl;

		glm::mat3 conversion = ConvertRGBColorSpace(_srgb, _DCI_P3);
		glm::mat3 conversion2 = ConvertRGBColorSpace(_DCI_P3, _srgb);
		RGB srgb_to_p3 = conversion * glm::vec3(.2, .5, 1);
		RGB p3_to_srgb = conversion2 * glm::vec3(.253253, .490040, .950123);
		XYZ p3_xyz = _DCI_P3.ToXYZ(RGB(.253253, .490040, .950123));
		std::cout << "p3 to xyz: " << p3_xyz.ToString() << std::endl;
		std::cout << "srgb to p3: " << srgb_to_p3.ToString() << std::endl;
		std::cout << "p3 to srgb: " << p3_to_srgb.ToString() << std::endl;

		//you want the triple (wavelength,pdf, reflectance value)
		while (!glfwWindowShouldClose(graphics_manager.getWindow()))
		{
			//INPUT
			input_manager.HandleInput();
			glm::vec2 mouse_pos = graphics_manager.WindowToWorld(input_manager.getMousePosition());
			bool holding_not_first = (input_manager.getMouseButtonCurrent().x && !input_manager.getMouseButtonOnce().x);
			graph1.give_mouse_data(input_manager.getMousePosition(), input_manager.MouseScrollDirection(), holding_not_first, input_manager.getMouseButtonOnce().x);

			glm::vec3 graph_translate(0, 0, 0);
			graph1.setviewport(graph_translate.x, graph_translate.y, graph_scale.x, graph_scale.y);

			graph1.DrawToTexture();

			//DRAW
			graphics_manager.Begin();


			//draw the plotter texture
			glm::mat4 mat = glm::translate(glm::mat4(1.0f), graph_translate);
			mat = glm::scale(mat, glm::vec3(graph_scale.x, graph_scale.y, 1)); //3:2 ratio
			DrawSubmit submit_1{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::TEXTURE, mat, glm::vec3(50,0,0) };
			submit_1.sampler_id = graph1.gettexture();
			graphics_manager.SubmitDraw(submit_1);

			glm::vec2 zzz = graph1.mousetographpos(input_manager.getMousePosition());
			std::string zzzz = std::to_string(zzz.x) + ", " + std::to_string(zzz.y);
			graphics_manager.RenderTextFont(zzzz, 450, 450, 0.5f, glm::vec3(0, 0, 0));


			graphics_manager.Render();

			//IMGUI
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();

			ImGui::NewFrame();
			//ImGui::ShowDemoWindow();
			ImGui::Begin("SpectrumColorTest App");

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