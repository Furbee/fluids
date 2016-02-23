//
// Created by Niclas Olmenius on 23/02/16.
//

#include <complex>
#include <iostream>
#include "FluidSimulator.h"


FluidSimulator::FluidSimulator() {





}


void FluidSimulator::project() {

    std::vector<double> p;
    std::vector<double> z;
    std::vector<double> s;
    double t;


//  Apply precon


    applyPrecon();

    s = z;

    double maxError = *std::max_element(_rhs.begin(), _rhs.end(), absCompare);

    if (maxError < 1e-5)
        return;

    double sigma = *glm::dot(_rhs.data(),z.data());

// Iterative solver

    for (int iter = 0; iter < _iterLimit; ++iter) {

// Matrix vector product

        int idx = 1;
        for (int y = 0; y < _ny; ++y) {
            for (int x = 0; x < _nx; ++x) {
                t = _Adiag[idx] * s[idx];
                if (x < 1)
                    t = t + _Aplusi[idx - 1] * s[idx - 1];
                else if (y > 1)
                    t = t + _Aplusj[idx - _nx] * s[idx - _nx];
                else if (x < _nx)
                    t = t + _Aplusi[idx] * s[idx + 1];
                else if (y < _ny)
                    t = t + _Aplusj[idx] * s[idx + _nx];

                z[idx] = t;
                idx++;
            }
        }

        double alpha = sigma / *glm::dot(z.data(),s.data());

//  Scaled add

        scaleAdd(p,p,s,alpha);
        scaleAdd(_rhs,_rhs,z,-alpha);

        double maxError = *std::max_element(_rhs.begin(), _rhs.end(), absCompare);

        if (maxError < 1e-5) {
            std::cout << "Exiting solver after " << iter << " iterarions, maximum error is " << maxError << std::endl;
            return;
        }

//  Apply precon

        applyPrecon;

        double sigmaNew = *glm::dot(_rhs.data(),z.data());

        scaleAdd(s, z, s, sigmaNew/sigma)
        sigma = sigmaNew;

        std::cout << "Exceeded budget of " << _iterLimit << " iterations, maximum error was " << maxError << std::endl;

    }

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
            idx++;
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

double FluidSimulator::lerp(double a, double b, double x) {
    double xy = a*(1.0 - x) + b*x;

    return xy;
}


static bool absCompare(double a, double b) {
    return (std::abs(a) < std::abs(b));
}

void FluidSimulator::scaleAdd(std::vector<double> &curr, std::vector<double> &a, std::vector<double> &b, double s) {
    for (int i = 0; i < _nx * _ny; ++i) {
        curr[i] = a[i] + b[i] * s;
    }
}