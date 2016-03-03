//
// Created by Niclas Olmenius on 23/02/16.
//

#include "Renderer.h"
#include "Shader.h"
#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>
#include <OpenGL/gl.h>

using namespace glm;


Renderer::Renderer(unsigned int width, unsigned int height) : IMAGE_WIDTH(width), IMAGE_HEIGHT(height) {
    frameBuffer = 0;

}


void Renderer::error_callback(int error, const char *description) {
    fputs(description, stderr);
}

void Renderer::key_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}


void Renderer::init(unsigned char* image) {

    _image = image;

    glfwSetErrorCallback(error_callback);

    //initiate GLFW
    glfwInit();


    if (!glfwInit()) {
        std::cout << "Initialization failed" << std::endl;
        exit(EXIT_FAILURE);
    }
    //Set options for GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

    //Create window
    window = glfwCreateWindow(WIDTH, HEIGHT, "Fluid", nullptr, nullptr);
    glfwMakeContextCurrent(window);
    if (window == nullptr) {

        std::cout << "I am a failure, for even a simple window I can not open." << std::endl;
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwSetKeyCallback(window, key_callback);

    glewExperimental = GL_TRUE;
    glewInit();

    glViewport(0, 0, WIDTH, HEIGHT);


    shader = new Shader("../shaders/vertShader.vert", "../shaders/fragShader.frag");



//    glGenVertexArrays(1, &VAO);
//    glGenBuffers(1, &VBO);

//    glBindVertexArray(VAO);

//    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    GLuint frameBuffer = 0;
    glGenFramebuffers(1, &frameBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

    glGenTextures(1, &renderedTexture);
    // "Bind" the newly created texture : all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, renderedTexture);
    // Filtering
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    // Give an empty image to OpenGL ( the last "0" )
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, IMAGE_WIDTH, IMAGE_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, _image);

    glBindTexture(GL_TEXTURE_2D,0);


//    glBindBuffer(GL_ARRAY_BUFFER, 0);

//    glBindVertexArray(0);


}

void Renderer::render() {

    glfwPollEvents();

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);



    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, renderedTexture);
    glUniform1i(glGetUniformLocation(shader->Program, "ourTexture1"), 0);
    shader->Use();



    glfwSwapBuffers(window);
    //glfwWaitEvents();


}

void Renderer::cleanup() {

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);


    glfwDestroyWindow(window);

    glfwTerminate();
    exit(EXIT_SUCCESS);


}
