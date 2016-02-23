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


    double alpha = (_DENSITYSOOT - _DENSITYAIR)/_DENSITYAIR;

    for (int y = 0, idx = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++) {

            double buoyancy = _dt * _GRAVITY * (alpha * _d[idx] - (_T(idx) - TAMB) / TAMB);
            _v[idx] += buoyancy * 0.5;
            _v[idx + _nx] += buoyancy * 0.5;

        }
    }


}

void FluidSimulator::advect() {

    int idx = 0;
    double ix, iy;
    struct Point x0y0;

    for( int y = 0 ; y < _ny ; y++ ){
        for( int x = 0 ; y < _nx ; x++ ){

            //offset
            ix = x + 0.5;
            iy = y + 0.5;

            x0y0 = //rungeKutta(ix, iy, dt, _u, _v, _dxy, _nx, _ny);
             // = //rungeKutta(ix, iy, dt, _u, _v, _dxy, _nx, _ny);

            _dn[idx] = cerp2( x0y0.x, x0y0.y, _nx, _ny, 0.5, 0.5, _d  );

            _dn[idx] = std::max(0, _dn[idx] * exp( -KDISS * _dt) );

            _Tn[idx] = cerp2( x0y0.x, x0y0.y, _nx, _ny, 0.5, 0.5, _T );

            idx++;

        }

    }

    idx = 0;

    for( int y = 0; y < _ny ; y++ ){
        for(int x = 0 ; x <= _nx ; x++){

            //offset
            ix = x + 0.5;
            iy = y + 0.5;

            x0y0.x = //rungeKutta(ix, iy, dt, _u, _v, _dxy, _nx, _ny);
            x0y0.y = //rungeKutta(ix, iy, dt, _u, _v, _dxy, _nx, _ny);

            _un[idx] = cerp2( x0y0.x, x0y0.y, (_nx + 1), _ny, 0.0, 0.5, _u );
            idx++

        }

    }

    idx = 0;
    for( int y = 0 ; y <= _ny ; y++ ){
        for( int x = 0 ; x < _nx ; x++ ) {

            //offset
            ix = x + 0.5;
            iy = y + 0.5;

            x0y0.x = //rungeKutta(ix, iy, dt, _u, _v, _dxy, _nx, _ny);
            x0y0.y = //rungeKutta(ix, iy, dt, _u, _v, _dxy, _nx, _ny);

            _vn[idx] = cerp2( x0y0.x, x0y0.y, _nx, (_ny+1), 0.5, 0.0, _v);
            idx++;
        }

    }
    //update u & v
    _d = _dn;
    _T = _Tn;
    _u = _un;
    _v = _vn;
}
