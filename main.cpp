
#include "FluidSimulator.h"
#include "Renderer.h"


int main(void) {

    const unsigned int texWidth = 128, texHeight = 128;

    FluidSimulator fluidSimulator;
    Renderer renderer(texWidth,texHeight);
    renderer.init(fluidSimulator.getDensityImage());

    double time = 0.0;

    while (time < 4.0) {

        fluidSimulator.update();
        renderer.render();

        time += fluidSimulator.getTimestep();

    }

    renderer.cleanup();



}
