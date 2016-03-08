//
// Created by Niclas Olmenius on 23/02/16.
//

#include <complex>
#include <iostream>
#include <numeric>
#include <iomanip>
#include <algorithm>

#include "FluidSimulator.h"


FluidSimulator::FluidSimulator(unsigned int width, unsigned int height) {

    _nx = width;
    _ny = height;
    _nz = 32;

    length = _nx * _ny;

    _dxy = 1.0 / std::min(_nx, _ny);
    _dt = 0.0025;
//    _dt = 1.0 / 60.0;
    _umax = 0.0;

    std::srand(std::time(0));

    _Adiag = std::vector<double>(length, 0);
    _Aplusi = std::vector<double>(length, 0);
    _Aplusj = std::vector<double>(length, 0);
    _Aplusk = std::vector<double>(length, 0);
    _rhs = std::vector<double>(length, 0);
    _pressure = std::vector<double>(length, 0);
    _s = std::vector<double>(length, 0);
    _z = std::vector<double>(length, 0);
    _u = std::vector<double>((_nx + 1) * _ny, 0);
    _un = std::vector<double>((_nx + 1) * _ny, 0);
    _v = std::vector<double>(_nx * (_ny + 1), 0);
    _vn = std::vector<double>(_nx * (_ny + 1), 0);
//    std::vector<double> _w
//    std::vector<double> _wn;

    _d = std::vector<double>(length, 0);
    _dn = std::vector<double>(length, 0);
    _T = std::vector<double>(length, 0);
    _Tn = std::vector<double>(length, 0);
    _precon = std::vector<double>(length, 0);


    _image = new unsigned char[_nx * _ny * 4];


}

FluidSimulator::~FluidSimulator() {

    delete _image;

}

void FluidSimulator::project() {

    double t;

    std::fill(_pressure.begin(), _pressure.end(), 0.0);
    std::fill(_z.begin(), _z.end(), 0.0);

//  Apply precon

    applyPrecon();

    _s = _z;

    double maxError = std::abs(*std::max_element(_rhs.begin(), _rhs.end(), absCompare));

    if (maxError < 1e-5)
        return;

    double sigma = std::inner_product(_rhs.begin(), _rhs.end(), _z.begin(), 0.0);

// Iterative solver

    for (int iter = 0; iter < _ITERLIMIT; ++iter) {

// Matrix vector product

        for (int y = 0, idx = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++, idx++) {
                t = _Adiag[idx] * _s[idx];

                if (x > 0) {
                    t = t + _Aplusi[idx - 1] * _s[idx - 1];
                }
                if (y > 0) {
                    t = t + _Aplusj[idx - _nx] * _s[idx - _nx];
                }
                if (x < _nx - 1) {
                    t = t + _Aplusi[idx] * _s[idx + 1];
                }
                if (y < _ny - 1) {
                    t = t + _Aplusj[idx] * _s[idx + _nx];
                }

                _z[idx] = t;
            }
        }

        double alpha = sigma / std::inner_product(_z.begin(), _z.end(), _s.begin(), 0.0);

//  Scaled add

        scaleAdd(_pressure, _pressure, _s, alpha);
        scaleAdd(_rhs, _rhs, _z, -alpha);

        maxError = std::abs(*std::max_element(_rhs.begin(), _rhs.end(), absCompare));

        if (maxError < 1e-5) {
            std::cout << "Exiting solver after " << iter << " iterarions, maximum error is " << maxError << std::endl;
            return;
        }

//  Apply precon

        applyPrecon();

        double sigmaNew = std::inner_product(_rhs.begin(), _rhs.end(), _z.begin(), 0.0);

        scaleAdd(_s, _z, _s, sigmaNew / sigma);
        sigma = sigmaNew;

    }

    std::cout << "Exceeded budget of " << _ITERLIMIT << " iterations, maximum error was " << maxError << std::endl;


}

void FluidSimulator::addInFlow(double x0, double y0, double x1, double y1, int w, int h, double ox, double oy,
                               double dxy,
                               double value,
                               std::vector<double> &src) {

    int ix0 = (int) (x0 / dxy - ox);
    int iy0 = (int) (y0 / dxy - oy);
    int ix1 = (int) (x1 / dxy - ox);
    int iy1 = (int) (y1 / dxy - oy);


    for (int y = std::max(iy0, 0); y < std::min(iy1, h); y++) {
        for (int x = std::max(ix0, 0); x < std::min(ix1, h); x++) {
            double l = hypot((2.0 * (x + 0.5) * dxy - (x0 + x1)) / (x1 - x0),
                             (2.0 * (y + 0.5) * dxy - (y0 + y1)) / (y1 - y0));

            double temp = std::min(std::abs(l), 1.0);
            double vi = (1.0 - temp * temp * (3.0 - 2.0 * temp)) * value;

            int idx = getIdx(x, y, w);
            if (std::abs(src[idx]) < std::abs(vi)) {
                src[idx] = vi;
            }

        }
    }


}

void FluidSimulator::update() {

//
//    double maxV = *std::max_element(_v.begin(), _v.end(), absCompare);
//    double maxU = *std::max_element(_u.begin(), _u.end(), absCompare);
//
//    double umax = std::max(maxU, maxV + sqrt(5 * _dxy * std::abs(_GRAVITY)));
//    _dt = umax > 0.00005 ? _dxy/umax : 0.0025;


//    std::cout << "addInFlow" << std::endl;



    int random_variable = 75 + std::rand() % 30;

    double inflow_x = 0.45;
    double inflow_x_width = 0.10;
    double inflow_y = 0.80;
    double inflow_y_height = 0.03;


    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
              _nx, _ny, 0.5, 0.5, _dxy, 1.0, _d);

    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
              _nx, _ny, 0.5, 0.5, _dxy, TAMB + 300.0, _T);

    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
              _nx + 1, _ny, 0.0, 0.5, _dxy, 5.0 * cos(random_variable * PI / 180), _u);

    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
              _nx, _ny + 1, 0.5, 0.0, _dxy, -1.0, _v);

//    inflow_x = 0.15;
//    inflow_y = 0.80;
//
//    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
//              _nx, _ny, 0.5, 0.5, _dxy, 1.0, _d);
//
//    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
//              _nx, _ny, 0.5, 0.5, _dxy, TAMB + 300.0, _T);
//
//    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
//              _nx + 1, _ny, 0.0, 0.5, _dxy, 5.0 * cos(random_variable * PI / 180), _u);
//
//    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
//              _nx, _ny + 1, 0.5, 0.0, _dxy, -1.0, _v);
//
//    inflow_x = 0.75;
//    inflow_y = 0.80;
//
//    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
//              _nx, _ny, 0.5, 0.5, _dxy, 1.0, _d);
//
//    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
//              _nx, _ny, 0.5, 0.5, _dxy, TAMB + 300.0, _T);
//
//    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
//              _nx + 1, _ny, 0.0, 0.5, _dxy, 5.0 * cos(random_variable * PI / 180), _u);
//
//    addInFlow(inflow_x, inflow_y, inflow_x + inflow_x_width, inflow_y + inflow_y_height,
//              _nx, _ny + 1, 0.5, 0.0, _dxy, -1.0, _v);



//    std::cout << "buildRhs" << std::endl;
    buildRhs();
//    std::cout << "PressureMatrix" << std::endl;
    buildPressureMatrix();
//    std::cout << "buildPrecon" << std::endl;
    buildPrecon();
//    std::cout << "project" << std::endl;
    project();
//    std::cout << "applyPressure" << std::endl;
    applyPressure();
//    std::cout << "advect" << std::endl;
    advect();
//    std::cout << "Buoyancy" << std::endl;
    applyBuoyancy();


}

void FluidSimulator::applyPressure() {

    // Apply pressure

    double scale = _dt / (_RHO * _dxy);
    int uvidx = 0;

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++, idx++) {
            uvidx = getIdx(x, y, _nx + 1);
            _u[uvidx] -= scale * _pressure[idx];
            _u[uvidx + 1] += scale * _pressure[idx];
            uvidx = getIdx(x, y, _nx);
            _v[uvidx] -= scale * _pressure[idx];
            _v[uvidx + _nx] += scale * _pressure[idx];
        }
    }

    // Update boundaries

    unsigned int idx = 0;

    for (int y = 0; y < _ny; y++) {
        idx = getIdx(0, y, _nx + 1);
        _u[idx] = 0.0;
        _u[idx + 1] = 0.0;
        idx = getIdx(_nx - 1, y, _nx + 1);
        _u[idx] = 0.0;
        _u[idx + 1] = 0.0;

        idx = getIdx(0, y, _nx);
        _T[idx] = TAMB;
        idx = getIdx(_nx - 1, y, _nx);
        _T[idx] = TAMB;

    }


    for (int x = 0; x < _nx; x++) {
        idx = getIdx(x, 0, _nx);
        _v[idx] = 0.0;
        _v[idx + 1] = 0.0;
        idx = getIdx(x, _ny - 1, _nx);
        _v[idx] = 0.0;
        _v[idx + _nx] = 0.0;

        idx = getIdx(x, 0, _nx);
        _T[idx] = TAMB;
        idx = getIdx(x, _ny - 1, _nx);
        _T[idx] = TAMB;
    }

//
//    for (int y = 0, ; y < _ny; y++, idx++) {
//        idx = getIdx(0, y, _nx);
//        _T[idx] = TAMB;
//        idx = getIdx(_nx - 1, y, _nx);
//        _T[idx] = TAMB;
//    }

//
//    for (int x = 0, idx = 0; x < _nx; x++, idx++) {
//        idx = getIdx(x, 0, _nx);
//        _T[idx] = TAMB;
//        idx = getIdx(x, _ny - 1, _nx);
//        _T[idx] = TAMB;
//    }
}

void FluidSimulator::buildRhs() {

    double scale = 1.0 / _dxy;

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

            if (e < _sigma_mic * _Adiag[idx]) {
                e = _Adiag[idx];
            }

            _precon[idx] = 1.0 / sqrt(e);

        }

    }
}

void FluidSimulator::applyPrecon() {

    for (int y = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++) {
            unsigned int idx = getIdx(x, y, _nx);

            double t = _rhs[idx];

            if (x > 0) {
                t -= _Aplusi[idx - 1] * _precon[idx - 1] * _z[idx - 1];
            }

            if (y > 0) {
                t -= _Aplusj[idx - _nx] * _precon[idx - _nx] * _z[idx - _nx];
            }

            _z[idx] = t * _precon[idx];

        }
    }

    for (int y = _ny - 1; y >= 0; y--) {
        for (int x = _nx - 1; x >= 0; x--) {

            unsigned int idx = getIdx(x, y, _nx);

            double t = _z[idx];

            if (x < _nx - 1) {
                t -= _Aplusi[idx] * _precon[idx] * _z[idx + 1];
            }
            if (y < _ny - 1) {
                t -= _Aplusj[idx] * _precon[idx] * _z[idx + _nx];
            }

            _z[idx] = t * _precon[idx];

        }
    }
}

void FluidSimulator::applyBuoyancy() {


    double alpha = (_DENSITYSOOT - _DENSITYAIR) / _DENSITYAIR;

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++, idx++) {

            double buoyancy = _dt * _GRAVITY * (alpha * _d[idx] - (_T[idx] - TAMB) / TAMB);

            _v[idx] += buoyancy * 0.5;
            _v[idx + _nx] += buoyancy * 0.5;

        }
    }

    for (int x = 0; x < _nx; x++) {
        int idx = getIdx(x, 0, _nx);
        _v[idx] = 0.0;
        _v[idx + 1] = 0.0;
        idx = getIdx(x, _ny - 1, _nx);
//        std::cout << idx << std::endl;
        _v[idx] = 0.0;
        _v[idx + _nx] = 0.0;
    }


}


void FluidSimulator::advect() {
    int index = 0;

    for (int idY = 0; idY < _ny; idY++) {
        for (int idX = 0; idX < _nx; idX++) {

            //offset
            double x = idX + 0.5;
            double y = idY + 0.5;

            rungeKutta3(x, y);

            _dn[index] = cerp2(x, y, _nx, _ny, 0.5, 0.5, _d);

            _dn[index] = std::max(0.0, _dn[index] * exp(-KDISS * _dt));

            _Tn[index] = cerp2(x, y, _nx, _ny, 0.5, 0.5, _T);
            index++;
        }
    }

    index = 0;

    for (int idY = 0; idY < _ny; idY++) {
        for (int idX = 0; idX <= _nx; idX++) {

            //offset
            double x = idX + 0.0;
            double y = idY + 0.5;

            rungeKutta3(x, y);

            _un[index] = cerp2(x, y, (_nx + 1), _ny, 0.0, 0.5, _u);
            index++;

        }
    }

    index = 0;

    for (int idY = 0; idY <= _ny; idY++) {
        for (int idX = 0; idX < _nx; idX++) {

            //offset
            double x = idX + 0.5;
            double y = idY + 0.0;

            rungeKutta3(x, y);

            _vn[index] = cerp2(x, y, _nx, (_ny + 1), 0.5, 0.0, _v);
            index++;

        }
    }

    //update the vectors
    _d = _dn;
    _T = _Tn;
    _u = _un;
    _v = _vn;
}


void FluidSimulator::rungeKutta3(double &x, double &y) {
    double earlyU = lerp2(x, y, 0.0, 0.5, _nx + 1, _ny, _u) / _dxy;
    double earlyV = lerp2(x, y, 0.5, 0.0, _nx, _ny + 1, _v) / _dxy;

    double mX = x - (0.5 * _dt * earlyU);
    double mY = y - (0.5 * _dt * earlyV);

    double mU = lerp2(mX, mY, 0.0, 0.5, _nx + 1, _ny, _u) / _dxy;
    double mV = lerp2(mX, mY, 0.5, 0.0, _nx, _ny + 1, _v) / _dxy;

    double lateX = x - (0.75 * _dt * mU);
    double lateY = y - (0.75 * _dt * mV);

    double lateU = lerp2(lateX, lateY, 0.0, 0.5, _nx + 1, _ny, _u);
    double lateV = lerp2(lateX, lateY, 0.5, 0.0, _nx, _ny + 1, _v);

    x -= _dt * ((2.0 / 9.0) * earlyU + (3.0 / 9.0) * mU + (4.0 / 9.0) * lateU);
    y -= _dt * (((2.0 / 9.0) * earlyV + (3.0 / 9.0) * mV + (4.0 / 9.0) * lateV));

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

    double xy = a * (1.0 - x) + b * x;

    return xy;
}

double FluidSimulator::cerp(double a, double b, double c, double d, double x) {

    double xsq = x * x;
    double xcu = xsq * x;
    double minV = std::min(a, std::min(b, std::min(c, d)));
    double maxV = std::max(a, std::max(b, std::max(c, d)));

    double t =
            a * (0.0 - 0.5 * x + 1.0 * xsq - 0.5 * xcu) +
            b * (1.0 + 0.0 * x - 2.5 * xsq + 1.5 * xcu) +
            c * (0.0 + 0.5 * x + 2.0 * xsq - 1.5 * xcu) +
            d * (0.0 + 0.0 * x - 0.5 * xsq + 0.5 * xcu);

    return std::min(std::max(t, minV), maxV);

}


void FluidSimulator::scaleAdd(std::vector<double> &curr, std::vector<double> &a, std::vector<double> &b, double s) {
    for (int i = 0; i < _nx * _ny; i++) {
        curr[i] = a[i] + b[i] * s;
    }
}

void FluidSimulator::buildPressureMatrix() {

    // Set up matrix entities for the pressure equations
    double scale = _dt / (_RHO * _dxy * _dxy);
    std::fill(_Adiag.begin(), _Adiag.end(), 0.0); // Coeffcient matrix for pressure equations
    std::fill(_Aplusi.begin(), _Aplusi.end(), 0.0);
    std::fill(_Aplusj.begin(), _Aplusj.end(), 0.0);

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++, idx++) {

            if (x < _nx - 1) {
                _Adiag[idx] += scale;
                _Adiag[idx + 1] += scale;
                _Aplusi[idx] = (-scale);
            } else {
                _Aplusi[idx] = 0.0;
            }

            if (y < _ny - 1) {
                _Adiag[idx] += scale;
                _Adiag[idx + _nx] += scale;
                _Aplusj[idx] = (-scale);
            } else {
                _Aplusj[idx] = 0.0;
            }

        }
    }


}

double FluidSimulator::lerp2(double x, double y, double ox, double oy, int w, int h,
                             const std::vector<double> &quantity) {
    x = std::min(std::max(x - ox, 0.0), w - 1.001);
    y = std::min(std::max(y - oy, 0.0), h - 1.001);

    int ix = (int) (x);
    int iy = (int) (y);
    x = x - ix;
    y = y - iy;


    double x00 = quantity[getIdx(ix, iy, w)];
    double x10 = quantity[getIdx(ix + 1, iy, w)];
    double x01 = quantity[getIdx(ix, iy + 1, w)];
    double x11 = quantity[getIdx(ix + 1, iy + 1, w)];

    return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);


}

double FluidSimulator::getTimestep() {
    return _dt;
}


void FluidSimulator::updateImage() {

    for (int i = 0; i < _nx * _ny; i++) {
        unsigned char shade = (unsigned char) ((1.0 - _d[i]) * 255.0);
        shade = std::max(std::min(shade, (unsigned char) 255), (unsigned char) 0);

        _image[i * 4 + 0] = shade;
        _image[i * 4 + 1] = shade;
        _image[i * 4 + 2] = shade;
        _image[i * 4 + 3] = 0xFF;

    }
}

unsigned char *FluidSimulator::getImagePtr() {
    return _image;
}


