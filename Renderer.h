//
// Created by Niclas Olmenius on 23/02/16.
//


#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include "Shader.h"


class Renderer {

public:
    Renderer();

    static void draw();


private:

    //Functions

    static void error_callback(int error, const char *description);

    static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);

    static void render(GLFWwindow *window);

};


