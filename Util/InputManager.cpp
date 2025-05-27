#include "../pch.h"
#include "InputManager.h"

InputManager::InputManager()
{
	mouse_pos = glm::vec2(0, 0);
	mouse_button_current = glm::bvec2(0, 0);
	mouse_button_once = glm::bvec2(0, 0);
	button_once_ready = glm::bvec2(1, 1);
	button_once_frame = glm::bvec2(0, 0);
	mouse_scroll = glm::vec2(0, 0);
	scrolled = false;

	mouse_middle_button = false;
	mouse_middle_once = false;
	mouse_middle = false;

	for (int i = 0; i < 500; i++) {
		Keys[i] = false;
	}
}

void InputManager::HandleInput()
{
	glfwPollEvents();

	if (scrolled)
		scrolled = false;
	else
		mouse_scroll = glm::vec2(0, 0);

	if (mouse_middle_once)
		mouse_middle_once = false;
	else
		mouse_middle_button = false;

	//mouse buttons
	if (mouse_button_once.x)
	{
		if (button_once_frame.x == false)
		{
			button_once_frame.x = true;
		}
		else
		{
			mouse_button_once.x = false;
			button_once_frame.x = false;
		}
	}

	if (mouse_button_once.y)
	{
		if (button_once_frame.y == false)
		{
			button_once_frame.y = true;
		}
		else
		{
			mouse_button_once.y = false;
			button_once_frame.y = false;
		}
	}
}

void InputManager::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS) {
		SetKeys(key, true);
		SetKeysOnce(key, true);
	}
	else if (action == GLFW_RELEASE) {
		SetKeys(key, false);
		SetKeysOnce(key, false);
	}
}

void InputManager::cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	mouse_pos = glm::vec2(xpos, ypos);
}

bool InputManager::GetKeyPress(int i) const
{
	return Keys[i];
}

bool InputManager::GetKeyPressOnce(int i)
{
	bool val = KeysOnce[i];
	KeysOnce[i] = false;
	return val;
}

void InputManager::mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_3)
	{
		if (action == GLFW_PRESS)
		{
			mouse_middle_button = true;
			mouse_middle_once = true;

			mouse_middle = true;
		}
		else if (action == GLFW_RELEASE) {
			mouse_middle = false;
		}
	}

	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		if (action == GLFW_PRESS) {
			mouse_button_current.x = true;
			if (button_once_ready.x == true)
			{
				mouse_button_once.x = true;
				button_once_ready.x == false;
			}
		}
		else if (action == GLFW_RELEASE) {
			mouse_button_current.x = false;
			button_once_ready.x == true;
		}
	}

	if (button == GLFW_MOUSE_BUTTON_RIGHT) {
		if (action == GLFW_PRESS) {
			mouse_button_current.y = true;
			if (button_once_ready.y == true)
			{
				mouse_button_once.y = true;
				button_once_ready.y == false;
			}
		}
		else if (action == GLFW_RELEASE) {
			mouse_button_current.y = false;
			button_once_ready.y == false;
		}
	}
}

void InputManager::scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	mouse_scroll = glm::vec2(xoffset, yoffset);
	scrolled = true;
}