//
// Created by Niclas Olmenius on 23/02/16.
//

#include "FluidSimulator.h"


FluidSimulator::FluidSimulator() {





}


void FluidSimulator::project() {

}

void FluidSimulator::addInFlow(int x0, int y0, int x1, int y1, int w, int h, int ox, int oy, double dxy, double value,
                               std::unique_ptr<std::vector<double>> src) {

}

void FluidSimulator::update(double dt) {

}

void FluidSimulator::buildRhs() {

    double scale = 1/_dxy;

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++, idx++) {
            _rhs[idx] = -scale * ((_u[getIdx(x+1,y,_nx+1)] - _u[getIdx(x,y,_nx+1)])
                        + (_v[getIdx(x,y+1,_ny)] - _v[getIdx(x,y,_ny)]));

        }

    }

}
