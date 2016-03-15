#include "FluidSimulator.h"
#include "Renderer.h"
#include <iomanip>
#include "external/lodepng/lodepng.h"


int main(void) {

    // Size of simulation and the generated texture
    const unsigned int texWidth = 128, texHeight = 128;

    FluidSimulator fluidSimulator(texWidth, texHeight);
    Renderer renderer(texWidth, texHeight);
    renderer.init(fluidSimulator.getImagePtr()); // Init the renderer

    double time = 0.0;

    unsigned int iterations = 0;
    double iter_time = 0.0;

    while (time < 20.0 && !glfwWindowShouldClose(renderer.window)) {

        // Run the simulation for 0.016 seconds between every exported frame
        while (iter_time < 0.016) {
            fluidSimulator.update();
            time += fluidSimulator.getTimestep();
            fluidSimulator.updateImage();
            renderer.render();
            iter_time += fluidSimulator.getTimestep();
        }

        iter_time = 0.0;

        // Generate png image
        std::ostringstream ostring_stream;
        ostring_stream << "images/Frame" << std::setfill('0') << std::setw(5) << iterations++ << ".png";
        lodepng::encode(ostring_stream.str(), fluidSimulator.getImagePtr(), texWidth, texHeight);


    }

    renderer.cleanup();


}
