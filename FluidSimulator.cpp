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

void FluidSimulator::applyPressure() {

    // Apply pressure

    double scale = _dt / (_rho * _dxy);
    int uvidx = 0;
    int idx = 0;

    for (int y = 0; y < _ny; ++y) {
        for (int x = 0; x < _nx; ++x) {
            uvidx = getIdx(x, y, _nx + 1);
            _u[uvidx] = _u[uvidx] - scale * _u[idx];
            _u[uvidx + 1] = _u[uvidx + 1] + scale * _u[idx];
            uvidx = getIdx(x, y, _nx);
            _v[uvidx] = _v[uvidx] - scale * _u[idx];
            _v[uvidx + _nx] = _v[uvidx + _nx] + scale * _u[idx];
            idx = idx++;
        }
    }

    // Update boundaries

    for (int y = 0; y < _ny; ++y) {
        idx = getIdx(1, y, _nx + 1);
        _u[idx] = 0.0;
        _u[idx + 1] = 0.0;
        idx = getIdx(_nx, y, _nx + 1);
        _u[idx] = 0.0;
        _u[idx + 1] = 0.0;
    }


    for (int x = 0; x < _nx; ++x) {
        idx = getIdx(x, 1, _nx);
        _v[idx] = 0.0;
        _v[idx + 1] = 0.0;
        idx = getIdx(x, _ny, _nx);
        _v[idx] = 0.0;
        _v[idx + _nx] = 0.0;
    }


    for (int y = 0; y < _ny; ++y) {
        idx = getIdx(1, y, _nx);
        _T[idx] = tAmb;
        idx = getIdx(_nx, y, _nx);
        _T[idx] = tAmb;
    }


    for (int x = 0; x < _nx; ++x) {
        idx = getIdx(x, 1, _nx);
        _T[idx] = tAmb;
        idx = getIdx(x, _ny, _nx);
        _T[idx] = tAmb;
    }
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

void FluidSimulator::buildPrecon() {

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++, idx++) {

            double e = _Adiag[idx];

            if (x > 0) {
                double px = _Aplusi[idx - 1]*_precon[idx - 1];
                double py = _Aplusj[idx - 1]*_precon[idx - 1];

                e = e - (px*px + _tau_mic*px*py);
            }
            if (y > 0) {
                double px = _Aplusi[idx - _nx]*_precon[idx - _nx];
                double py = _Aplusj[idx - _nx]*_precon[idx - _nx];
                e = e - (py*py + _tau_mic*px*py);
            }

            if (e < _sigma_mic*_Adiag[idx])
                e = _Adiag[idx];

            _precon[idx] = 1.0/sqrt(e);

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


void FluidSimulator::applyBuoyancy() {


    double alpha = (_densitySoot - _densityAir)/_densityAir;

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++) {

            double buoyancy = _dt * _gravity * (alpha * _d[idx] - (_T[idx] - tAmb) / tAmb);
            _v[idx] += buoyancy * 0.5;
            _v[idx + _nx] += buoyancy * 0.5;

        }
    }


}

double FluidSimulator::cerp2(double x, double y, int w, int h, double ox, double oy, std::vector<double> &quantity) {

    x = std::min(std::max(x - ox, 0.0), w - 1.001);
    y = std::min(std::max(y - oy, 0.0), h - 1.001);
    int ix = static_cast<int>(x);
    int iy = static_cast<int>(y);
    x = x - ix;
    y = y - iy;

    int x0 = std::max(ix - 1, 0);
    int x1 = ix;
    int x2 = ix + 1;
    int x3 = std::min(ix + 2, w - 1);

    int y0 = std::max(iy - 1, 0);
    int y1 = iy;
    int y2 = iy + 1;
    int y3 = std::min(iy + 2, h - 1);

    double q0 = cerp(quantity[getIdx(x0, y0, w)], quantity[getIdx(x1, y0, w)],
    quantity[getIdx(x2, y0, w)], quantity[getIdx(x3, y0, w)], x);
    double q1 = cerp(quantity[getIdx(x0, y1, w)], quantity[getIdx(x1, y1, w)],
    quantity[getIdx(x2, y1, w)], quantity[getIdx(x3, y1, w)], x);
    double q2 = cerp(quantity[getIdx(x0, y2, w)], quantity[getIdx(x1, y2, w)],
    quantity[getIdx(x2, y2, w)], quantity[getIdx(x3, y2, w)], x);
    double q3 = cerp(quantity[getIdx(x0, y3, w)], quantity[getIdx(x1, y3, w)],
    quantity[getIdx(x2, y3, w)], quantity[getIdx(x3, y3, w)], x);

    return cerp(q0, q1, q2, q3, y);

}

double FluidSimulator::lerp(double a, double b, double x) {

    double xy = a*(1.0 - x) + b*x;

    return xy;
}

double FluidSimulator::cerp(double a, double b, double c, double d, double x) {

    double xsq = x*x;
    double xcu = xsq*x;
    double minV = std::min(a, std::min(b, std::min(c, d)));
    double maxV = std::max(a, std::max(b, std::max(c, d)));

    double t =
    a*(0.0 - 0.5*x + 1.0*xsq - 0.5*xcu) +
    b*(1.0 + 0.0*x - 2.5*xsq + 1.5*xcu) +
    c*(0.0 + 0.5*x + 2.0*xsq - 1.5*xcu) +
    d*(0.0 + 0.0*x - 0.5*xsq + 0.5*xcu);

    return std::min(std::max(t, minV), maxV);

}
