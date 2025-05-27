#pragma once
class InputManager
{
public:
	InputManager();

	//call every frame
	void HandleInput();

	void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
	void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
	void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
	void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

	glm::vec2 getMousePosition() { return mouse_pos; }
	glm::bvec2 getMouseButtonCurrent() { return mouse_button_current; }
	glm::bvec2 getMouseButtonOnce() { return mouse_button_once; }
	glm::vec2 getMouseScroll() { return mouse_scroll; }
	bool getMouseMiddleDown() { return mouse_middle_button; }
	bool getMouseMiddleCurrent() { return mouse_middle; }
	int MouseScrollDirection()
	{
		return mouse_scroll.y;
	}

	void SetKeys(int i, bool val) { Keys[i] = val; }
	void SetKeysOnce(int i, bool val) { KeysOnce[i] = val; }

	bool* GetKeys() { return &Keys[0]; }
	bool GetKeyPress(int i) const;
	bool GetKeyPressOnce(int i);
private:
	//mouse
	glm::vec2 mouse_pos;
	glm::bvec2 mouse_button_current;
	glm::bvec2 mouse_button_once;
	bool mouse_middle;
	bool mouse_middle_button;
	bool mouse_middle_once;
	glm::vec2 mouse_scroll;
	bool scrolled;

	glm::bvec2 button_once_ready;
	glm::bvec2 button_once_frame;

	//keys
	bool Keys[500];
	bool KeysOnce[500];
};

