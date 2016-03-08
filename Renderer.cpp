//
// Created by Niclas Olmenius on 23/02/16.
//

#include "Renderer.h"

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


void Renderer::init(unsigned char *image) {

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

    // fix for HiDPI displays
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);

    shader = new Shader("../shaders/vertShader.vert", "../shaders/fragShader.frag");


    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);


    // data for a fullscreen quad (this time with texture coords)
    GLfloat vertexData[] = {
            //  X     Y     Z           U     V
            1.0f, 1.0f, 0.0f, 1.0f, 1.0f, // vertex 0
            -1.0f, 1.0f, 0.0f, 0.0f, 1.0f, // vertex 1
            1.0f, -1.0f, 0.0f, 1.0f, 0.0f, // vertex 2
            -1.0f, -1.0f, 0.0f, 0.0f, 0.0f, // vertex 3
    }; // 4 vertices with 5 components (floats) each

    // fill with data
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 4 * 5, vertexData, GL_STATIC_DRAW);


    // set up generic attrib pointers
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (char *) 0 + 0 * sizeof(GLfloat));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), (char *) 0 + 3 * sizeof(GLfloat));


    // generate and bind the index buffer object
    glGenBuffers(1, &IBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);

    GLuint indexData[] = {
            0, 1, 2, // first triangle
            2, 1, 3, // second triangle
    };

    // fill with data
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * 2 * 3, indexData, GL_STATIC_DRAW);

    // "unbind" vao
    glBindVertexArray(0);


    // generate a texture
    glGenTextures(1, &renderedTexture);
    // "Bind" the newly created texture : all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, renderedTexture);


    // Filtering
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // Give an empty image to OpenGL ( the last "0" )
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, IMAGE_WIDTH, IMAGE_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);




}

void Renderer::render() {




    glClear(GL_COLOR_BUFFER_BIT);


    shader->Use();

    // activate and bind our texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, renderedTexture);
    // upload new image data to the gpu
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, IMAGE_WIDTH, IMAGE_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, _image);

    // set the uniform of our current shader program
    glUniform1i(glGetUniformLocation(shader->Program, "ourTexture1"), 0);

    // bind vertex array
    glBindVertexArray(VAO);

    // draw
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);


    // check for errors
    GLenum error = glGetError();
    if (error != GL_NO_ERROR) {
        std::cerr << error << std::endl;

    }


    glfwSwapBuffers(window);
    //glfwWaitEvents();
    glfwPollEvents();


}

void Renderer::cleanup() {

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);


    glfwDestroyWindow(window);

    glfwTerminate();
    exit(EXIT_SUCCESS);


}
