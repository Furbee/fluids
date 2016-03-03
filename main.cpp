#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "FluidSimulator.h"
#include "Renderer.h"



#define GLSL(src) #src


int main(void) {

    FluidSimulator fluidSimulator;

    double time = 0.0;

    while(time < 8.0) {

        fluidSimulator.update();

        time += fluidSimulator.getTimestep();

    }




    //Renderer::draw();

}
