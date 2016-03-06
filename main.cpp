#include "FluidSimulator.h"
#include "Renderer.h"


int main(void) {

    const unsigned int texWidth = 128, texHeight = 128;

    FluidSimulator fluidSimulator(texWidth, texHeight);
    Renderer renderer(texWidth, texHeight);
    renderer.init(fluidSimulator.getImagePtr());

    double time = 0.0;

    while (time < 20.0) {

        for (int i = 0; i < 4; i++) {
            fluidSimulator.update();
            time += fluidSimulator.getTimestep();
        }
        fluidSimulator.updateImage();
        renderer.render();


    }

    renderer.cleanup();


}
