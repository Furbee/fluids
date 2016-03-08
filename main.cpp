#include "FluidSimulator.h"
#include "Renderer.h"
#include <iomanip>
#include "external/lodepng/lodepng.h"


int main(void) {

    const unsigned int texWidth = 128, texHeight = 128;

    FluidSimulator fluidSimulator(texWidth, texHeight);
    Renderer renderer(texWidth, texHeight);
    renderer.init(fluidSimulator.getImagePtr());

    double time = 0.0;

    unsigned int iterations = 0;

    while (time < 20.0) {

        for (int i = 0; i < 4; i++) {
            fluidSimulator.update();
            time += fluidSimulator.getTimestep();
            fluidSimulator.updateImage();
            renderer.render();
        }

        std::ostringstream ostring_stream;
        ostring_stream << "images/Frame" << std::setfill('0') << std::setw(5) << iterations++ << ".png";
        lodepng::encode(ostring_stream.str(), fluidSimulator.getImagePtr(), texWidth, texHeight);


    }

    renderer.cleanup();


}
