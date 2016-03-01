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


void Renderer::error_callback(int error, const char *description) {
    fputs(description, stderr);
}

void Renderer::key_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

void Renderer::render(GLFWwindow *window) {

}

void Renderer::draw() {

    /*******************************************
     *              SIZE OF WINDOW
     ******************************************/
    const GLuint HEIGHT = 480, WIDTH = 640;

    /*******************************************
     *              DECLARATIONS
     ******************************************/

    GLFWwindow *window;

    Shader Shaders("shaders/vertShader.vert", "shaders/fragShader.frag");


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

    // Vertex corner points
    GLfloat vertices[] = {
            -0.5f, -0.5f, 0.0f, //LEFT
            0.5f, 0.5f, 0.0f,   //RIGHT
            0.0f, 0.5f, 0.0f    //TOP
    };
    GLuint VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    Shaders.Use();
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid *) 0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindVertexArray(0);

    /***********************************************
     *                  RUNTIME
     **********************************************/

    while (!glfwWindowShouldClose(window)) {
        render(window);

        glfwPollEvents();

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glBindVertexArray(0);

        glfwSwapBuffers(window);
        //glfwWaitEvents();
    }
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);


    glfwDestroyWindow(window);

    glfwTerminate();
    exit(EXIT_SUCCESS);

}

