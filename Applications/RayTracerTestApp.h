#pragma once
#include "../Graphics/GraphicsManager.h"
#include "../RayTracer/Film.h"
#include "../ThirdParty/pbrv4/samplers.h"
#include "../ThirdParty/pbrv4/filters.h"
//#include "../ThirdParty/pbrv4/"

#include "../RayTracer/Shapes.h"
#include "../RayTracer/Cameras.h"
#include "../RayTracer/AssetManager.h"
#include "../RayTracer/Octtree_Model.h"

#include <thread>
#include <condition_variable>
#include <atomic>


class RayTracerTestApp
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		RayTracerTestApp* handler = reinterpret_cast<RayTracerTestApp*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<RayTracerTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<RayTracerTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<RayTracerTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<RayTracerTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons;
		if (graphics_manager.InitWindow(1850, 1850 / 1.5f, "Ray Tracer Test App", icons) == EXIT_FAILURE)
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
		/* setup raytracing pipeline, loading data and creating everything, then do the whole pipeline, output to image
		*/

		//Models
		MeshCache::LoadMeshFromFile("Game_Data/models/teapot.fbx");
		MeshCache::LoadMeshFromFile("Game_Data/models/monkeyhead.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/stanford-bunny.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/stanford-dragon.obj");

		//Textures
		int bitDepth = -1;
		GLFWimage images[1];
		stbi_set_flip_vertically_on_load(true);
		images[0].pixels = stbi_load("Game_Data/uvmap.png", &images[0].width, &images[0].height, &bitDepth, 3);
		stbi_set_flip_vertically_on_load(false);

		//SHAPES
		glm::mat4 m1 = glm::translate(glm::mat4(1.0f), glm::vec3(0, -40, 800)) 
			         * glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0, 1, 0))
			         * glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(1, 0, 0))
			         * glm::scale(glm::mat4(1.0f), glm::vec3(15, 15, 15)); //15 for dragon
		Cylinder cylinder_1("cylinder1", m1, 100, -50, 50, 250.0);
		//Sphere sphere_1("sphere1", m1, 100, -100, 100, 360.0f);
		//TriangleSimple triangle_1("triangle1", m1, glm::vec3(-500, 0, -250), glm::vec3(500, 0, -250), glm::vec3(0, 0, 250));
		//Disk disk_1("disk1", m1, 0, 50, 150, 250.0f);
		//Triangle triangle_2("triangle2", m1, "Game_Data/models/monkeyhead.obj", 0, 0);

		TriModel model1("model1", m1, "Game_Data/models/stanford-dragon.obj", true, false, Triangle::vertex_available());
		model1.ComputeBackFace(glm::vec3(0, 0, 1), true);
		Octtree_Model oct_model1(model1);
		oct_model1.CreateOcttree();
		//oct_model1.PrintInfo();

		int shape_id = 0;
		std::vector<Shape*> Shapes = { &cylinder_1 };

		//
		int film_X = 500;
		int film_Y = 500;
		int film_resolution_X = 500;
		int film_resolution_Y = 500;

		size_t data_size = film_resolution_X * film_resolution_Y * 3;//x*y*rgb
		std::vector<unsigned char> buffer(data_size);
		for (int i = 0; i < buffer.size(); i += 3)
		{
			buffer[i] = 255;
			buffer[i + 1] = 0;
			buffer[i + 2] = 0;
		}
		//to start use the basic pinhole camera, and just create 1 ray per pixel and do a 1 intersection test, output to opengl texture to draw a rect
		GLuint film_texture_id;
		glGenTextures(1, &film_texture_id);
		glBindTexture(GL_TEXTURE_2D, film_texture_id);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, film_resolution_X, film_resolution_Y, 0, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE_EXT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE_EXT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

		/*
		std::vector<glm::vec3> pixel_values(data_size / 3);
		for (int i = 0; i < pixel_values.size(); i++)
		{
			if (i == 0)
				pixel_values[i] = glm::vec3(255, 0, 0);
			else
				pixel_values[i] = glm::vec3(0, 0, 0);
		}
		*/

		pbrt::Spectra::Init();
		pbrt::RGBToSpectrumTable::Init();
		pbrt::RGBColorSpace::Init();

		//filters
		glm::vec2 pixel_size(film_X / (float)film_resolution_X, film_Y / (float)film_resolution_Y);
		std::cout << "pixel size: " << pixel_size.x << " " << pixel_size.y << std::endl;
		pbrt::BoxFilter box_filter(pixel_size / glm::vec2(2, 2));
		pbrt::GaussianFilter gauss_filter(pixel_size / glm::vec2(2,2));
		pbrt::TriangleFilter tri_filter(pixel_size / glm::vec2(2, 2));
		//pixel sensor
		pbrt::PixelSensor sensor_xyz(pbrt::RGBColorSpace::sRGB, pbrt::GetNamedSpectrum("stdillum-D65"), 1.0f / pbrt::CIE_Y_integral);

		//doesnt work
		pbrt::PixelSensor sensor_canon(pbrt::GetNamedSpectrum("canon_eos_100d_r"), pbrt::GetNamedSpectrum("canon_eos_100d_g"), 
			pbrt::GetNamedSpectrum("canon_eos_100d_b"), pbrt::RGBColorSpace::sRGB, pbrt::GetNamedSpectrum("stdillum-A"), 1.0f / pbrt::CIE_Y_integral);
		

		Film film_1;
		film_1.film_dim = glm::ivec2(film_X, film_Y);
		film_1.image_res = glm::ivec2(film_resolution_X, film_resolution_Y);
		film_1.filter = &tri_filter;
		film_1.pixel_sensor = &sensor_xyz;
		film_1.pixels.resize(film_resolution_X * film_resolution_Y);
		
		//camera
		glm::vec3 cam_pos(0, 0, 0);
		float cam_pitch = 90.0f; //up,down starting at (0,1,0)
		float cam_yaw = 90.0f;
		float speed = 1.0f;
		float fov = 45.0f;
		//PinholeCamera pin_cam(0, glm::vec3(film_X, film_Y, 10), glm::vec3(0, 0, 0), glm::vec3(0, 0, 1), glm::vec3(1, 0, 0), glm::vec3(0, 1, 0), glm::vec2(film_resolution_X, film_resolution_Y));

		//OrthographicCamera ortho_cam(0, 1, film_X, film_Y, glm::vec3(0, 0, 0), glm::vec3(0, 0, 1), glm::vec3(1, 0, 0), glm::vec3(0, 1, 0),
			//glm::vec2(film_resolution_X, film_resolution_Y));

		float lens_radius = 50;
		float focal_distance = 800.0f;
		PerspectiveCamera persp_cam(1.0f, 1000.0f, (float)film_1.film_dim.x, (float)film_1.film_dim.y, fov, glm::vec3(0, 0, 0), glm::vec3(0, 0, 1),
			glm::vec3(1, 0, 0), glm::vec3(0, 1, 0), glm::vec2(film_1.image_res.x, film_1.image_res.y), lens_radius, focal_distance); //(10,200)

		int cam_id = 0;
		std::vector<CameraBase*> cameras = { &persp_cam };
		CameraBase* main_cam = &persp_cam;

		//uniform random sampler
		int samples_per_pixel_pos = 10;
		int seed = 0;
		pbrt::IndependentSampler sampler_independant(samples_per_pixel_pos, seed);
		pbrt::StratifiedSampler strat_sampler(samples_per_pixel_pos, samples_per_pixel_pos, seed);

		pbrt::Sampler* main_sampler = &strat_sampler;
		std::cout << "starting\n";

		pbrt::Spectrum* illumA = pbrt::GetNamedSpectrum("stdillum-A");
		pbrt::Spectrum* illumD = pbrt::GetNamedSpectrum("stdillum-D65"); 
		pbrt::Spectrum* illumF = pbrt::GetNamedSpectrum("stdillum-F1");
		
		//THREADS
		bool spawn_threads = true;
		std::vector<std::thread> thread_pool;
		std::vector<pbrt::Sampler*> thread_samplers;
		std::condition_variable cv;
		std::atomic<int> jobs_done = 0;

		bool geometry_test = false;
		int pixel_index = 0;
		int index_limit = samples_per_pixel_pos;
		bool paused = false;
		float colors[3] = { 0.5,0.5,0.5 };
		while (!glfwWindowShouldClose(graphics_manager.getWindow()))
		{
			//INPUT
			input_manager.HandleInput();
			glm::vec2 mouse_pos = graphics_manager.WindowToWorld(input_manager.getMousePosition());

			persp_cam.SetlensRadius(lens_radius);
			persp_cam.SetfocalDistance(focal_distance);


			auto Li = [&](Ray ray, pbrt::SampledWavelengths lambdas) -> pbrt::SampledSpectrum
			{
				std::optional<LocalSurfaceInfo> surfaceoptional = surfaceoptional = oct_model1.Traverse(ray); //Shapes[shape_id]->Intersect(ray);
				if (surfaceoptional.has_value())
				{
					//world coordinate
					LocalSurfaceInfo surfaceinfo = surfaceoptional.value();
					glm::vec3 world_n = surfaceinfo.n;
					glm::vec3 world_p = surfaceinfo.hitp;

					pbrt::SampledSpectrum radiance = pbrt::SampledSpectrum(0);
					if (false)
					{
						//uv texel- (0,0) is the bottomleft pixel and we scan right-up
						int x_pixel = std::floor(surfaceinfo.u * images[0].width);
						int y_pixel = std::floor(surfaceinfo.v * images[0].height);
						int index = (y_pixel * images[0].width) + (x_pixel % images[0].width);
						glm::vec3 texel(images[0].pixels[index * 3 + 0], images[0].pixels[index * 3 + 1], images[0].pixels[index * 3 + 2]);
						texel = texel / glm::vec3(255, 255, 255);
						//convert rgb to spectral
						pbrt::RGBAlbedoSpectrum rgb_to_sprectrum(*pbrt::RGBColorSpace::sRGB, pbrt::RGB(texel));
						pbrt::SampledSpectrum texture_spectral = rgb_to_sprectrum.Sample(lambdas);

						radiance += texture_spectral;
					}
					else
					{
						//I think the spectral values have to be normalized or else, if its reflectance. These represent reflectance we can scale later
						pbrt::RGBIlluminantSpectrum light_spectrum = pbrt::RGBIlluminantSpectrum(*pbrt::RGBColorSpace::sRGB, pbrt::RGB(1,1,1));
						pbrt::SampledSpectrum light_spectral = light_spectrum.Sample(lambdas);
						
						//ambient - spectral of illuminants, scale correctly
						//dynamic_cast<pbrt::PiecewiseLinearSpectrum*>(illumF)->Scale(1 / illumF->MaxValue());
						pbrt::SampledSpectrum ambient_spectral = 0.3f * illumF->Sample(lambdas); //not normalized

						//object spectral reflectance
						pbrt::RGBAlbedoSpectrum mat_spectrum(*pbrt::RGBColorSpace::sRGB, pbrt::RGB(colors[0], colors[1], colors[2]));
						pbrt::SampledSpectrum mat_spectral = mat_spectrum.Sample(lambdas);

						//cosine
						float cosine_1 = glm::clamp(glm::dot(world_n, glm::vec3(0, 0, -1)), 0.0f, 1.0f);
						float cosine_2 = glm::clamp(glm::dot(world_n, glm::vec3(0, 1, 0)), 0.0f, 1.0f);
						float cosine_3 = glm::clamp(glm::dot(world_n, glm::vec3(1, -1, 1)), 0.0f, 1.0f);
						//glm::vec3 light_dir = glm::rotate(glm::mat4(1.0f), glm::radians(lightangle), glm::vec3(0, 1, 0)) * glm::vec4(0, 0, 1, 0);
						float light_1_cos = glm::clamp(glm::dot(world_n, glm::vec3(0, 0, -1)), 0.0f, 1.0f);
						

						radiance += ambient_spectral;
						radiance += light_1_cos * (light_spectral * mat_spectral);
					}

					if (false)
					{
						dynamic_cast<pbrt::PiecewiseLinearSpectrum*>(illumF)->Scale(1 / illumF->MaxValue());
						pbrt::SampledSpectrum ambient_spectral = 1.0f * illumF->Sample(lambdas); //not normalized

						//spotlight
						float spotlight_cos = glm::clamp(glm::dot(glm::normalize(cameras[cam_id]->getlookdirection()), glm::normalize(world_p - cam_pos)), 0.0f, 1.0f);
						if (spotlight_cos >= .9999)
							radiance += std::powf(spotlight_cos, 1) * ambient_spectral;
					}
					//pbrt::SampledSpectrum lambda_weights = lambdas.PDF();
					return radiance;// *lambdas.PDF();
				}		

				return pbrt::SampledSpectrum(0);
			};

			//index should just change random variables so it gets a new ray direction and such
			auto evaluate_pixel = [&](int pixel_id, int index, pbrt::Sampler* sampl) -> void
			{
				int x_pix = pixel_id % film_1.image_res.x;
				int y_pix = film_1.image_res.y - std::floor(pixel_id / (float)film_1.image_res.x);
				glm::ivec2 pixel(x_pix, y_pix);

				if (geometry_test)
				{
					Ray ray = main_cam->generateRay(pixel, sampl);
					pbrt::SampledSpectrum L = Li(ray, pbrt::SampledWavelengths());
					if (L != pbrt::SampledSpectrum(0))
					{
						film_1.pixels[pixel_id].rgbsum = glm::vec3(1, 0, 0);
						film_1.pixels[pixel_id].weightsum = 1;
					}
					return;
				}

				sampl->StartPixelSample(pixel, index, 0);

				//sample wavelengths
				//these are the wavelengths were going to pick the spectral values from
				pbrt::SampledWavelengths lambdas = pbrt::SampledWavelengths::SampleVisible(sampl->Get1D());

				//sample film position/time/camera stuff
				
				//calculate imagespace position (pixel + .5 + offset)

				//pixel offset is first done randomly with the main samplers, then its importance sampled based on the filter
				glm::vec2 uniform_pixel_offset = sampl->GetPixel2D();
				pbrt::FilterSample fs = film_1.filter->Sample(uniform_pixel_offset);
				//std::cout << fs.p.x << " " << fs.p.y << std::endl;
				glm::vec2 pixel_sampled_pos = glm::fvec2(pixel) + glm::vec2(.5f, .5f) + fs.p; //image space
			
				//might sample direction(if pinhole then its 1 direction and not random, if lens then it will uniform sample lens over 0,1 so no pdf)
				//pretty much our direction we pick we dont need a weight
				Ray ray = main_cam->generateRay(pixel_sampled_pos, sampl); //convert to film position

				//get radiance(sampled spectrum) from ray
				pbrt::SampledSpectrum L = Li(ray, lambdas); //scale by cam_direction weight 


				//convert it to cameraRGB-->outputRGB and add sample to pixel value
				pbrt::RGB cam_RGB = film_1.pixel_sensor->ToSensorRGB(L, lambdas);

				cam_RGB.r = glm::clamp(cam_RGB.r, 0.0f, 1.0f);
				cam_RGB.g = glm::clamp(cam_RGB.g, 0.0f, 1.0f);
				cam_RGB.b = glm::clamp(cam_RGB.b, 0.0f, 1.0f);
	
				film_1.pixels[pixel_id].rgbsum += fs.weight * cam_RGB.getglm(); //weight by position
				film_1.pixels[pixel_id].weightsum += fs.weight;




				//position of the pixel is weighted at the end of the rgb
				//the ray direction is weighted at the radiance of the ray returned
				//the wavelength is weighted anywhere we need to do those color space integral
			};

			if (!paused)
			{				
				auto ThreadFunction = [&](int begin_index, int end_index, pbrt::Sampler* sampl)
				{
					//lock doesn't do anything
					std::mutex mut;
					std::unique_lock<std::mutex> lck(mut);
	
					int offset = std::hash<std::thread::id>{}(std::this_thread::get_id()) % 255;
					float theoffset = offset / 255.0f;
					//std::cout << theoffset << std::endl;
					while (true)
					{
						cv.wait(lck);
						for (int j = begin_index; j < end_index; j++)
						{
							//film_1.pixels[j].rgbsum = glm::vec3(theoffset, theoffset, theoffset);
							//film_1.pixels[j].weightsum = 1.0f;
							evaluate_pixel(j, pixel_index, sampl);
						}
						jobs_done.fetch_add(1);
						//std::cout << jobs_done << std::endl;
					}
				};

				if (spawn_threads)
				{
					std::cout << "SPAWNING THREADS\n";
					int thread_count = std::max(std::thread::hardware_concurrency(), 1u);
					//thread_count = 1;
					int pixels_per_thread = film_1.pixels.size() / thread_count;
					thread_pool.resize(thread_count);
					thread_samplers.resize(thread_count);
					int beggining_id = 0;
					std::cout << "thread size: " << thread_count << " pixels_per_thread: " << pixels_per_thread << std::endl;
					for (int i = 0; i < thread_pool.size(); i++)
					{
						int end_id = beggining_id + pixels_per_thread;
						if (i == thread_pool.size() - 1)
							end_id = film_1.pixels.size();
						std::cout << "thread: " << i << " begggining: " << beggining_id << " end: " << end_id << std::endl;
						//pbrt::IndependentSampler samplerz_independant(samples_per_pixel_pos, seed);
						//pbrt::StratifiedSampler stratz_sampler(samples_per_pixel_pos, samples_per_pixel_pos, seed);
						pbrt::StratifiedSampler* stratz_sampler = new pbrt::StratifiedSampler(samples_per_pixel_pos, samples_per_pixel_pos, seed);
						thread_samplers[i] = stratz_sampler;
						thread_pool[i] = std::thread(ThreadFunction, beggining_id, end_id, stratz_sampler);
						beggining_id = end_id;
					}
					spawn_threads = false;
					std::this_thread::sleep_for(std::chrono::seconds(1));
				}
				//std::cout << "waking up\n";
				auto t1 = std::chrono::high_resolution_clock::now();
				cv.notify_all();
				//std::cout << "waiting\n";
				//jobs_done.wait(thread_pool.size());
				while (jobs_done.load() != thread_pool.size())
				{
					std::this_thread::sleep_for(std::chrono::milliseconds(1));
				}
				jobs_done = 0;
				auto t2 = std::chrono::high_resolution_clock::now();
				std::cout << "time: " << t2 - t1 << std::endl;
			//	std::cout << "finished\n";

				/*
				//pick pixel and index and gather pixel data
				for (int j = 0; j < film_1.pixels.size(); j++)
				{
					evaluate_pixel(j, pixel_index, main_sampler);
				}
				*/

				pixel_index++;
				if (pixel_index > index_limit)
					paused = true;

				//divide out all the pixel values to get final pixel result and convert to output rgb space rgb-->xyz->srgb
				for (int i = 0; i < film_1.pixels.size(); i++)
				{
					if (geometry_test)
					{
						float scale = 255.0f;
						//now write to gpu texture
						buffer[3 * i] = scale * film_1.pixels[i].rgbsum.r;
						buffer[3 * i + 1] = scale * film_1.pixels[i].rgbsum.g;
						buffer[3 * i + 2] = scale * film_1.pixels[i].rgbsum.b;
						continue;
					}

					pbrt::RGB sensor_rgb(film_1.pixels[i].rgbsum / film_1.pixels[i].weightsum);

					pbrt::XYZ xyz_val = film_1.pixel_sensor->XYZFromSensorRGB * sensor_rgb.getglm();
					pbrt::RGB output_rgb = pbrt::RGBColorSpace::sRGB->ToRGB(xyz_val);
					//output_rgb = sensor_rgb;

					output_rgb.r = glm::clamp(output_rgb.r, 0.0f, 1.0f);
					output_rgb.g = glm::clamp(output_rgb.g, 0.0f, 1.0f);
					output_rgb.b = glm::clamp(output_rgb.b, 0.0f, 1.0f);

					float scale = 255.0f;
					//now write to gpu texture
					buffer[3 * i] = scale * output_rgb.r;
					buffer[3 * i + 1] = scale * output_rgb.g;
					buffer[3 * i + 2] = scale * output_rgb.b;
				}

				glBindTexture(GL_TEXTURE_2D, film_texture_id);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, film_resolution_X, film_resolution_Y, 0, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());
			}

			//DRAW
			graphics_manager.Begin();


			glm::mat4 mat = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 0));
			mat = glm::scale(mat, glm::vec3(graphics_manager.getWidth(), graphics_manager.getHeight(), 1));
			DrawSubmit submit_1{ GEOMETRY_TYPE::RECTANGLE, SHADER_TYPE::TEXTURE, mat, glm::vec3(0,0,0) };
			submit_1.sampler_id = film_texture_id;
			graphics_manager.SubmitDraw(submit_1);

			graphics_manager.Render();

			
			//IMGUI
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();
			ImGui::NewFrame();
			//ImGui::ShowDemoWindow();
			ImGui::Begin("RayTracerTest app");

			ImGui::DragFloat("lens radius", &lens_radius);
			ImGui::DragFloat("focal distance", &focal_distance);

			if (ImGui::Button("restart"))
			{
				pixel_index = 0;
				for (int i = 0; i < film_1.pixels.size(); i++)
				{
					film_1.pixels[i].rgbsum = glm::vec3(0, 0, 0);
					film_1.pixels[i].weightsum = 0;
				}
				paused = false;
			}

			if (ImGui::Button("pause"))
			{
				paused = !paused;
			}

			if (ImGui::ColorPicker3("rgb to spectral", colors));
			colors[0] = glm::clamp(colors[0], 0.0f, 1.0f);
			colors[1] = glm::clamp(colors[1], 0.0f, 1.0f);
			colors[2] = glm::clamp(colors[2], 0.0f, 1.0f);

			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
				
			ImGui::End();
			ImGui::Render();
			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

			graphics_manager.End();

			//char a;
			//std::cin >> a;
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