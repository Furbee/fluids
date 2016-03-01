//
// Created by Oscar Nord on 01/03/16.
//
// Taken from learnopengl.com

#ifndef FLUID_SIMULATION_SHADER_H
#define FLUID_SIMULATION_SHADER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <GL/glew.h>

class Shader
{
public:
    GLuint Program;

    Shader(const GLchar* vertexPath, const GLchar* fragmentPath);
    // Uses the current shader
    void Use()
    {
        glUseProgram(this->Program);
    }
};

#endif //FLUID_SIMULATION_SHADER_H
