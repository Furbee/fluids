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

    /*
     * Initialize the GLFW context
     */
    void init(unsigned char *image);

    /*
     * Render the image using OpenGL
     */
    void render();

    /*
     * Cleanup the GLFW context
     */
    void cleanup();

    GLFWwindow *window; // Pointer to the GLFW window


private:

    /*
     * Constants
     */
    const GLuint HEIGHT = 480, WIDTH = 640;


    /*
     * Member variables
     */
    const GLuint IMAGE_WIDTH, IMAGE_HEIGHT;
    unsigned char *_image;
    GLuint VBO, VAO, IBO, renderedTexture;
    Shader *shader;


    /*
     * GLFW specific functions
     */
    static void error_callback(int error, const char *description);

    static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);


};


