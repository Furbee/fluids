//
// Created by Niclas Olmenius on 23/02/16.
//


#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <vector>
#include "Shader.h"


class Renderer {

public:
    Renderer(unsigned int width, unsigned int height);

    void init(unsigned char* image);
    void render();
    void cleanup();


private:

    // Member variables

    const GLuint HEIGHT = 480, WIDTH = 640;
    const GLuint IMAGE_HEIGHT, IMAGE_WIDTH;
    unsigned char* _image;

    GLuint VBO, VAO, IBO, frameBuffer, renderedTexture;
    Shader* shader;
    GLFWwindow *window;


    //Functions

    static void error_callback(int error, const char *description);

    static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);


};


