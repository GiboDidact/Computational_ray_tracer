#pragma once

#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <span>
#include <array>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <functional>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <string>
#include <numbers>
#include <memory>

#include "GL/glew.h"
#include "GLFW/glfw3.h" 

#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#include <glm\glm.hpp>
#include <glm\gtc\matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/matrix_operation.hpp>

#include "Util/HelperFunctions.h"
#include "Util/InputManager.h"
#include "Util/Timer.h"

//#include "Graphics/GraphicsManager.h"


constexpr float OneMinusEpsilon = 1 - std::numeric_limits<float>::min();
constexpr float Pi = 3.14159265358979323846;
constexpr float InvPi = 0.31830988618379067154;
constexpr float Inv2Pi = 0.15915494309189533577;
constexpr float Inv4Pi = 0.07957747154594766788;
constexpr float PiOver2 = 1.57079632679489661923;
constexpr float PiOver4 = 0.78539816339744830961;
constexpr float Sqrt2 = 1.41421356237309504880;