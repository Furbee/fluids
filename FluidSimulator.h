//
// Created by Niclas Olmenius on 23/02/16.
//


#pragma once

#include <vector>
#include <glm/glm.hpp>

class FluidSimulator {


public:
    FluidSimulator();


private:


    // Constants

    double _gravity = 9.82; // gravity
    double _rho = 0.1; // density (1e3 for water, 1.3 for air)
    double _densitySoot = 0.11;
    double _densityAir = 0.1;
    double tAmb = 273;
    double kDiss = 0.1;


    // Member variables

    int _nx;
    int _ny;
    int _nz;

    double _dxy;
    double _dt;
    double _umax;

    std::vector<double> _Adiag;
    std::vector<double> _Aplusi;
    std::vector<double> _Aplusj;
    std::vector<double> _Aplusk;
    std::vector<double> _rhs;
    std::vector<double> _pressure;
    std::vector<double> _u;
    std::vector<double> _z;
    std::vector<double> _un;
    std::vector<double> _v;
    std::vector<double> _vn;
    std::vector<double> _w;
    std::vector<double> _wn;

    std::vector<double> _d;
    std::vector<double> _dn;
    std::vector<double> _T;
    std::vector<double> _Tn;




    // MIC
    double _tau_mic = 0.97;
    double _sigma_mic = 0.25;
    std::vector<double> _precon;


    // Functions

    inline int getIdx(int i, int j, int width) {
        return i + j * width;
    }


    double lerp(double a, double b, double x );

    double lerp3(int x, int y, int z, float ox, float oy, float oz, int w, int h,
                 std::unique_ptr<std::vector<double>> quantity );

    double cerp( double a, double b, double c, double d, double x );

    double cerp2(double x, double y, int w, int h, double ox, double oy, std::vector<double> &quantity);


    void addInFlow( int x0, int y0, int x1, int y1, int w, int h, int ox, int oy, double dxy, double value,
                    std::unique_ptr<std::vector<double>> src);


    void applyBuoyancy();

    void buildRhs();

    void buildPrecon();

    void applyPrecon();

    void applyPressure();

    void advect();




    void update(double dt);


    void project();





};


