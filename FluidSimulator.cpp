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

    double scale = 1 / _dxy;

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++, idx++) {
            _rhs[idx] = -scale * ((_u[getIdx(x + 1, y, _nx + 1)] - _u[getIdx(x, y, _nx + 1)])
                                  + (_v[getIdx(x, y + 1, _ny)] - _v[getIdx(x, y, _ny)]));

        }

    }

}

void FluidSimulator::buildPrecon() {

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++, idx++) {

            double e = _Adiag[idx];

            if (x > 0) {
                double px = _Aplusi[idx - 1] * _precon[idx - 1];
                double py = _Aplusj[idx - 1] * _precon[idx - 1];

                e = e - (px * px + _tau_mic * px * py);
            }
            if (y > 0) {
                double px = _Aplusi[idx - _nx] * _precon[idx - _nx];
                double py = _Aplusj[idx - _nx] * _precon[idx - _nx];
                e = e - (py * py + _tau_mic * px * py);
            }

            if (e < _sigma_mic * _Adiag[idx])
                e = _Adiag[idx];

            _precon[idx] = 1.0 / sqrt(e);

        }
    }

}

void FluidSimulator::applyPrecon() {

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++, idx++) {
            double t = _rhs[idx];

            if (x > 1) {
                t = t - _Aplusi[idx - 1] * _precon[idx - 1] * _z[idx - 1];
            }

            if (y > 1) {
                t = t - _Aplusj[idx - _nx] * _precon[idx - _nx] * _z[idx - _nx];
            }

            _z[idx] = t * _precon[idx];

            idx = idx + 1;

        }
    }

    for (int y = _ny - 1, idx = _nx * _ny - 1; y >= 0; y--) {
        for (int x = _nx - 1; x >= 0; x--, idx--) {
            double t = _z[idx];

            if (x < _nx) {
                t = t - _Aplusi[idx] * _precon[idx] * _z[idx + 1];
            }
            if (y < _ny) {
                t = t - _Aplusj[idx] * _precon[idx] * _z[idx + _nx];
            }

            _z[idx] = t * _precon[idx];
        }
    }
}

