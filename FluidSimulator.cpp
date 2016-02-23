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

void FluidSimulator::applyPressure(double dt, double rho, int dxy, std::vector<double> &u, std::vector<double> &v,
                                   std::vector<double> &p, std::vector<double> &T) {

    // Apply pressure

    double scale = dt/(rho*dxy);
    int uvidx = 0;
    int idx = 0;

    for (int y = 0; y < _ny; ++y) {
        for (int x = 0; x < _nx; ++x) {
            uvidx = getIdx(x,y,_nx+1);
            u[uvidx] = u[uvidx] - scale * p[idx];
            u[uvidx+1] = u[uvidx+1] + scale * p[idx];
            uvidx = getIdx(x,y,_nx);
            v[uvidx] = v[uvidx] - scale * p[idx];
            v[uvidx + _nx] = v[uvidx + _nx] + scale * p[idx];
            idx = idx++;
        }
    }

    // Update boundaries

    for (int y=0; y < _ny; ++y) {
        idx = getIdx(1,y,_nx+1);
        u[idx] =  0.0;
        u[idx+1] = 0.0;
        idx = getIdx(_nx,y,_nx+1);
        u[idx] = 0.0;
        u[idx+1] = 0.0;
    }


    for (int x = 0; x < _nx; ++x) {
        idx = getIdx(x,1,_nx);
        v[idx] = 0.0;
        v[idx+1] = 0.0;
        idx = getIdx(x,_ny,_nx);
        v[idx] = 0.0;
        v[idx+_nx] = 0.0;
    }


    for (int y = 0; y < _ny; ++y) {
        idx = getIdx(1,y,_nx);
        T[idx] = tAmb;
        idx = getIdx(_nx,y,_nx);
        T[idx] = tAmb;
    }


    for (int x = 0; x < _nx; ++x) {
        idx = getIdx(x,1,_nx);
        T[idx] = tAmb;
        idx = getIdx(x,_ny,_nx);
        T[idx] = tAmb;
    }

}
