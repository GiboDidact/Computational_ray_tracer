#pragma once
#include "../../Graphics/GraphicsManager.h"
#include "../../RayTracer/Shapes.h"
#include "../../RayTracer/Cameras.h"
#include "../../RayTracer/AssetManager.h"
#include "../../RayTracer/Octtree_Model.h"
#include <thread>

class ShapeTestApp
{
	static void framebuffer_size_callback(GLFWwindow* window, int width, int height)
	{
		ShapeTestApp* handler = reinterpret_cast<ShapeTestApp*>(glfwGetWindowUserPointer(window));
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		handler->GetGraphicsManager().framebuffer_size_callback(window, w, h);
	}
	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		reinterpret_cast<ShapeTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().key_callback(window, key, scancode, action, mods);
	}
	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		reinterpret_cast<ShapeTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().cursor_position_callback(window, xpos, ypos);
	}
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		reinterpret_cast<ShapeTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().mouse_button_callback(window, button, action, mods);
	}
	static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		reinterpret_cast<ShapeTestApp*>(glfwGetWindowUserPointer(window))->GetInputManager().scroll_callback(window, xoffset, yoffset);
	}

public:
	void Start()
	{
		std::vector<std::string> icons;
		if (graphics_manager.InitWindow(1850, 1850 / 1.5f, "ShapeTest App", icons) == EXIT_FAILURE)
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
		//for each shape: test intersection, normals, uv, bounds, maybe sampling 
		//for each camera: test ray generations, maybe even visualize the rays

		//load models
		MeshCache::LoadMeshFromFile("Game_Data/models/cone.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/teapot.fbx");
		MeshCache::LoadMeshFromFile("Game_Data/models/backpack/backpack.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/cyber/0.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/monkeyhead.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/stanford-bunny.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/stanford-dragon.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/Stairs_2.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/crashb.obj");
		MeshCache::LoadMeshFromFile("Game_Data/models/Grunt.dae");

		int film_X = 500;
		int film_Y = 500;
		int film_resolution_X = 200;
		int film_resolution_Y = 200;

		size_t data_size = film_resolution_X * film_resolution_Y * 3;//x*y*rgb
		std::vector<unsigned char> buffer(data_size);
		for (int i = 0; i < buffer.size(); i+=3)
		{
			buffer[i] = 255;
			buffer[i+1] = 0;
			buffer[i+2] = 0;
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
		//no filtering texelFetch()
		//pixel buffer object texture streaming

		std::vector<glm::vec3> pixel_values(data_size / 3);
		for (int i = 0; i < pixel_values.size(); i++)
		{
			if(i == 0)
				pixel_values[i] = glm::vec3(255, 0, 0);
			else
				pixel_values[i] = glm::vec3(0, 0, 0);
		}

		//UV Texture
		int bitDepth = -1;
		GLFWimage images[1];
		stbi_set_flip_vertically_on_load(true);
		images[0].pixels = stbi_load("Game_Data/uvmap.png", &images[0].width, &images[0].height, &bitDepth, 3);
		stbi_set_flip_vertically_on_load(false);

		
		//SHAPES
		glm::mat4 m1 = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 100)) * glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0, 1, 0))
			* glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(1, 0, 0))
			* glm::scale(glm::mat4(1.0f), glm::vec3(250, 250, 250));
		Cylinder cylinder_1("cylinder1", m1, 100, -50, 50, 250.0);
		Sphere sphere_1("sphere1", m1, 100, -100, 100, 360.0f);
		TriangleSimple triangle_1("triangle1", m1, glm::vec3(-500, 0, -250), glm::vec3(500, 0, -250), glm::vec3(0, 0, 250));
		Disk disk_1("disk1", m1, 0, 50, 150, 250.0f);
		Triangle triangle_2("triangle2", m1, "Game_Data/models/monkeyhead.obj", 0, 0);

		TriModel model1("model1", m1, "Game_Data/models/stanford-bunny.obj", true, false, Triangle::vertex_available());
		model1.ComputeBackFace(glm::vec3(0, 0, 1), true);
		Octtree_Model oct_model1(model1);
		oct_model1.CreateOcttree();
		oct_model1.PrintInfo();
		//return;

		//return;
		int shape_id = 0;
		std::vector<Shape*> Shapes = { &sphere_1, &cylinder_1, &triangle_1, &triangle_2, &disk_1, &model1 };

		//CAMERA
		PinholeCamera pin_cam(0, glm::vec3(film_X, film_Y, 10), glm::vec3(0, 0, 0), glm::vec3(0, 0, 1), glm::vec3(1, 0, 0), glm::vec3(0, 1, 0), glm::vec2(film_resolution_X, film_resolution_Y));

		OrthographicCamera ortho_cam(0, 1, film_X, film_Y, glm::vec3(0, 0, 0), glm::vec3(0, 0, 1), glm::vec3(1, 0, 0), glm::vec3(0, 1, 0),
			                         glm::vec2(film_resolution_X, film_resolution_Y));
		
		PerspectiveCamera persp_cam(1.0f, 1000.0f, (float)film_X, (float)film_Y, 45.f, glm::vec3(0, 0, 0), glm::vec3(0, 0, 1), glm::vec3(1, 0, 0), glm::vec3(0, 1, 0),
			glm::vec2(film_resolution_X, film_resolution_Y));

		int cam_id = 0;
		std::vector<CameraBase*> cameras = { &ortho_cam, &persp_cam, &pin_cam };

		glm::vec3 cam_pos(0, 0, 0);
		float cam_pitch = 90.0f; //up,down starting at (0,1,0)
		float cam_yaw = 90.0f;
		float speed = 1.0f;

		bool bounding_boxes = false;
		float fov = 45.0f;
		bool usespotlight = false;
		bool useuv = false;
		float angle = 0.0f;
		float lightangle = 180.0f;
		float pinhole_z = 1.0f;
		float cosine_cutoff = .99f;
		bool print_error = false;
		bool pause = false;
		bool first_time = true;//somehow its throttled
		bool back_face_cull = true;
		bool kd_enable = false;
		int kd_node = 0;
		while (!glfwWindowShouldClose(graphics_manager.getWindow()))
		{
			//INPUT
			input_manager.HandleInput();
			glm::vec2 mouse_pos = graphics_manager.WindowToWorld(input_manager.getMousePosition());

			if (input_manager.GetKeyPressOnce(GLFW_KEY_P))
			{
				pause = !pause;
			}

			bool moved = false;
			//I can add it so I move with the camera orientation not just globally x,y,z
			if (input_manager.GetKeyPress(GLFW_KEY_W))
			{
				cam_pos += persp_cam.getlookdirection()*speed;
				moved = true;
			}
			if (input_manager.GetKeyPress(GLFW_KEY_S))
			{
				cam_pos -= persp_cam.getlookdirection() * speed;
				moved = true;
			}
			if (input_manager.GetKeyPress(GLFW_KEY_A))
			{
				cam_pos -= persp_cam.getrightdirection() * speed;
				moved = true;
			}
			if (input_manager.GetKeyPress(GLFW_KEY_D))
			{
				cam_pos += persp_cam.getrightdirection() * speed;
				moved = true;
			}
			if (input_manager.GetKeyPress(GLFW_KEY_SPACE))
			{
				cam_pos += persp_cam.getupdirection() * speed;
				moved = true;
			}
			if (input_manager.GetKeyPress(GLFW_KEY_LEFT_CONTROL))
			{
				cam_pos -= persp_cam.getupdirection() * speed;
				moved = true;
			}
			if (input_manager.GetKeyPress(GLFW_KEY_LEFT))
			{
				cam_yaw++;
				moved = true;
			}
			if (input_manager.GetKeyPress(GLFW_KEY_RIGHT))
			{
				cam_yaw--;
				moved = true;
			}
			if (input_manager.GetKeyPress(GLFW_KEY_UP))
			{
				cam_pitch--;
				moved = true;
			}
			if (input_manager.GetKeyPress(GLFW_KEY_DOWN))
			{
				cam_pitch++;
				moved = true;
			}

			if (moved == true)
			{
				for (int i = 0; i < cameras.size(); i++)
				{
					cameras[i]->setWorldPos(cam_pos);
					cameras[i]->setyawpitch(cam_yaw, cam_pitch);
				}
			}

			glm::mat4 m1 = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, 400)) * glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(0, 1, 0))
				* glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(1, 0, 0))
				* glm::scale(glm::mat4(1.0f), glm::vec3(100, 100, 100));
			//Shapes[shape_id]->SetRigidTransform(m1);
			
			if (back_face_cull && dynamic_cast<TriModel*>(Shapes[shape_id]))
			{
				dynamic_cast<TriModel*>(Shapes[shape_id])->ComputeBackFace(cameras[cam_id]->getlookdirection(), true);
			}

			
			//opengltexture starts bottomleft for first index, my film plane starts bottom left, 
			//ray trace and figure out every pixels (r,g,b) value

			if (!pause)
			{
				//spawn a bunch of threads and each on takes a i value and generates a pixel
				//or spread it out between the x threads, they write to the memory, maybe they have their own and we translate after
				auto evaluate_pixel = [&](int i) -> void
				{
					//for (int i = 0; i < pixel_values.size(); i++)
					//{
					//generate the ray for the pixel
					//opengl reads first pixel at bottomleft, our imagespace start at topleft so flip y
					int x_pix = i % film_resolution_X;
					int y_pix = film_resolution_Y - std::floor(i / (float)film_resolution_X);
					//std::cout << x_pix<<" "<<y_pix << std::endl;
					Ray ray = cameras[cam_id]->generateRay(glm::vec2(x_pix, y_pix));

					//intersect with a bounds
					if (bounding_boxes)
					{
						if (kd_enable)
						{
							bool intersectbounds = oct_model1.GetNode(kd_node).bounds.IntersectP(ray);
							if (intersectbounds)
							{
								pixel_values[i] = glm::vec3(255, 0, 0);
								//continue;
								return;
							}
						}
						else
						{
							bool intersectbounds = Shapes[shape_id]->Bounds().IntersectP(ray);
							if (intersectbounds)
							{
								pixel_values[i] = glm::vec3(255, 0, 0);
								//continue;
								return;
							}
						}
					}

					if (false)
					{
						std::optional<LocalSurfaceInfo> info = oct_model1.Traverse(ray);
					
						if (info.has_value())
						{
							pixel_values[i] = glm::vec3(255, 0, 0);
							return;
						}
						else
						{
							pixel_values[i] = glm::vec3(0, 0, 0);
							return;
						}
					}
					//intersect with an object
					//Timer t1;
					//t1.Begin();
					std::optional<LocalSurfaceInfo> surfaceoptional;
					if(kd_enable)
						surfaceoptional = oct_model1.Traverse(ray);
					else
						surfaceoptional = Shapes[shape_id]->Intersect(ray);
					//std::cout << t1.getTimeNano() << std::endl;
					if (surfaceoptional.has_value())
					{
						//pixel_values[i] = glm::vec3(255, 0, 0);
						//return;

						//world coordinate
						LocalSurfaceInfo surfaceinfo = surfaceoptional.value();
						glm::vec3 world_n = surfaceinfo.n;
						glm::vec3 world_p = surfaceinfo.hitp;


						pixel_values[i] = glm::vec3(0, 0, 0);
						if (useuv)
						{
							if (surfaceinfo.u > 1)
								std::cout << surfaceinfo.u << std::endl;
							//uv texel- (0,0) is the bottomleft pixel and we scan right-up
							int x_pixel = std::floor(surfaceinfo.u * images[0].width);
							int y_pixel = std::floor(surfaceinfo.v * images[0].height);
							int index = (y_pixel * images[0].width) + (x_pixel % images[0].width);
							glm::vec3 texel(images[0].pixels[index * 3 + 0], images[0].pixels[index * 3 + 1], images[0].pixels[index * 3 + 2]);

							pixel_values[i] += texel;
						}
						else
						{
							//ambient
							glm::vec3 ambient(10, 10, 10);

							//cosine
							float cosine_1 = glm::clamp(glm::dot(world_n, glm::vec3(0, 0, -1)), 0.0f, 1.0f);
							float cosine_2 = glm::clamp(glm::dot(world_n, glm::vec3(0, 1, 0)), 0.0f, 1.0f);
							float cosine_3 = glm::clamp(glm::dot(world_n, glm::vec3(1, -1, 1)), 0.0f, 1.0f);
							glm::vec3 light_dir = glm::rotate(glm::mat4(1.0f), glm::radians(lightangle), glm::vec3(0, 1, 0)) * glm::vec4(0, 0, 1, 0);
							float cosine_light = glm::clamp(glm::dot(world_n, light_dir), 0.0f, 1.0f);

							pixel_values[i] = ambient;
							pixel_values[i] += cosine_light * glm::vec3(100, 100, 100);
						}

						if (usespotlight)
						{
							//spotlight
							float spotlight_cos = glm::clamp(glm::dot(glm::normalize(cameras[cam_id]->getlookdirection()), glm::normalize(world_p - cam_pos)), 0.0f, 1.0f);
							if (spotlight_cos >= cosine_cutoff)
								pixel_values[i] += std::powf(spotlight_cos, 4) * glm::vec3(100, 100, 100);
						}


						pixel_values[i].x = glm::clamp(pixel_values[i].x, 0.0f, 255.0f);
						pixel_values[i].y = glm::clamp(pixel_values[i].y, 0.0f, 255.0f);
						pixel_values[i].z = glm::clamp(pixel_values[i].z, 0.0f, 255.0f);
					}
					else
					{
						pixel_values[i] = glm::vec3(0, 0, 0);
					}
				};
				for (int i = 0; i < pixel_values.size(); i++)
				{
					evaluate_pixel(i);
				}
				//spawn all the threads at the start and then push the indices onto the queue and then they grab them from the queue and run and they
				//then say they completed it and I wait until everythings completed to move on
				/*
				int thread_count = std::thread::hardware_concurrency();
				int pixels_per_thread = pixel_values.size() / 8;
				std::vector<std::thread> threads(thread_count);
				for (int i = 0; i < pixel_values.size(); i++)
				{
					int index = (i + 1) % thread_count;
					threads[index](evaluate_pixel, i);
				}
				*/
				//std::thread t1(evaluate_pixel, 0);

			}
			first_time = false;
			std::cout << "tris intersected: "<< Hitdata::triangle_intersect_count << std::endl;
			Hitdata::triangle_intersect_count = 0;
			//write to the texture
			for (int i = 0; i < pixel_values.size(); i ++)
			{
				buffer[3*i] = pixel_values[i].x;
				buffer[3 * i + 1] = pixel_values[i].y;
				buffer[3 * i + 2] = pixel_values[i].z;
			}
			glBindTexture(GL_TEXTURE_2D, film_texture_id);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, film_resolution_X, film_resolution_Y, 0, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());

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
			ImGui::Begin("ShapeTest App");

			ImGui::DragFloat("pinhole z", &pinhole_z);
			if (ImGui::DragFloat("fov", &fov))
			{
				persp_cam.ChangeFOV(fov);
			}
			ImGui::DragFloat("angle", &angle);
			ImGui::DragFloat("lightangle", &lightangle);
			ImGui::DragFloat("spotlight cutoff", &cosine_cutoff, .1, 0.0f, 1.0f);

			ImGui::SliderInt("Shape ID", &shape_id, 0, Shapes.size() - 1);
			ImGui::SliderInt("Cam ID", &cam_id, 0, cameras.size()-1);

			ImGui::Checkbox("draw uv", &useuv);
			ImGui::Checkbox("draw spotlight", &usespotlight);
			ImGui::Checkbox("bounding boxes", &bounding_boxes);
			if (ImGui::Checkbox("cull backfaces", &back_face_cull))
			{
				if (dynamic_cast<TriModel*>(Shapes[shape_id]))
				{
					dynamic_cast<TriModel*>(Shapes[shape_id])->EnableBackface(back_face_cull);
				}
			}
			ImGui::Checkbox("kd enable", &kd_enable);
			ImGui::InputInt("kd_node", &kd_node);

			if (ImGui::Button("error"))
			{
				print_error = !print_error;
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