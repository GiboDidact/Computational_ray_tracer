#pragma once
#include "fstream"

class ShaderCompiler
{
public:
	static std::string ReadFile(const std::string& fileLocation)
	{
		std::string content;
		std::ifstream fileStream(fileLocation.c_str(), std::ios::in);

		if (!fileStream.is_open()) {
			std::cout << "Failed to read" << fileLocation << "! File doesn't exist.\n";
			return "";
		}

		std::string line = "";
		while (!fileStream.eof())
		{
			std::getline(fileStream, line);
			content.append(line + "\n");
		}

		fileStream.close();
		return content;
	}

	static void printShaderlog(unsigned int shader)
	{
		int len;
		char* log;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
		if (len > 0)
		{
			log = (char*)malloc(len);
			glGetShaderInfoLog(shader, len, 0, log);
			std::cout << "Shader info log: " << std::string(log) << "\n";
			free(log);
		}
	}

	static void printProgramlog(unsigned int program)
	{
		int len;
		char* log;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
		if (len > 0)
		{
			log = (char*)malloc(len);
			glGetProgramInfoLog(program, len, 0, log);
			std::cout << "Program info log: " << std::string(log) << "\n";
			free(log);
		}
	}

	static bool CompileGLShader(const std::string& vertex, const std::string& frag, GLuint& programid)
	{
		bool success = true;

		std::string shadersource = ReadFile(vertex);
		std::string fragsource = ReadFile(frag);

		const char* vshadersource = shadersource.c_str();
		const char* fshadersource = fragsource.c_str();

		unsigned int vShader = glCreateShader(GL_VERTEX_SHADER);
		unsigned int fShader = glCreateShader(GL_FRAGMENT_SHADER);

		glShaderSource(vShader, 1, &vshadersource, NULL);
		glShaderSource(fShader, 1, &fshadersource, NULL);
		glCompileShader(vShader);
		glCompileShader(fShader);

		GLint compiled;
		glGetShaderiv(vShader, GL_COMPILE_STATUS, &compiled);
		if (compiled != 1)
		{
			std::cout << "vertex shader failed to compile!\n";
			success = false;
			printShaderlog(vShader);
		}
		glGetShaderiv(fShader, GL_COMPILE_STATUS, &compiled);
		if (compiled != 1)
		{
			std::cout << "fragment shader failed to compile!\n";
			success = false;
			printShaderlog(fShader);
		}

		programid = glCreateProgram();
		glAttachShader(programid, vShader);
		glAttachShader(programid, fShader);
		glLinkProgram(programid);

		glGetProgramiv(programid, GL_LINK_STATUS, &compiled);
		if (compiled != 1)
		{
			std::cout << "Program linking failed!\n";
			success = false;
			printProgramlog(programid);
		}

		return success;
	}
};