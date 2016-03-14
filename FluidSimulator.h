//
// Created by Niclas Olmenius on 23/02/16.
//


#pragma once

#include <vector>
#include <cmath>

class FluidSimulator {


public:
    FluidSimulator(unsigned int, unsigned int);


    virtual ~FluidSimulator();

    void update();

    double getTimestep();

    unsigned char *getImagePtr();

    void updateImage();


private:


    // Constants

    const double _GRAVITY = 9.82; // gravity
    const double _RHO = 0.1; // density (1e3 for water, 1.3 for air)
    const double _DENSITYSOOT = 0.25;
    const double _DENSITYAIR = 0.1;
    const double TAMB = 273;
    const double KDISS = 0.001;
    const int _ITERLIMIT = 600;
    const double PI = 3.14159265;
    const double _TOL = 10e-5;

    int pulse_time;

    // Member variables

    unsigned char *_image;

    unsigned int _nx;
    unsigned int _ny;
    int _nz;

    unsigned int length;

    double _dxy;
    double _dt;
    double _umax;


    std::vector<double> _Adiag;
    std::vector<double> _Aplusi;
    std::vector<double> _Aplusj;
    std::vector<double> _Aplusk;
    std::vector<double> _rhs;
    std::vector<double> _pressure;
    std::vector<double> _s;
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

    std::vector<double> _precon;
    //struct arrays for x,y,z
    struct Point {
        double x, y, z;
    };


    // MIC
    double _tau_mic = 0.97;
    double _sigma_mic = 0.25;


    // Functions

    inline unsigned int getIdx(int i, int j, int width) {
        return i + j * width;
    }


    double lerp(double a, double b, double x);


    double lerp2(double x, double y, double ox, double oy, int w, int h, const std::vector<double> &quantity);


    //double lerp3(int x, int y, int z, float ox, float oy, float oz, int w, int h,
    //             std::unique_ptr<std::vector<double>> quantity);

    double cerp(double a, double b, double c, double d, double x);

    double cerp2(double x, double y, int w, int h, double ox, double oy, std::vector<double> &quantity);


    void addInFlow();

    void addInFlow(double x0, double y0, double x1, double y1, int w, int h, double ox, double oy, double dxy,
                   double value,
                   std::vector<double> &src);


    void applyBuoyancy();

    void buildRhs();

    void buildPressureMatrix();

    void buildPrecon();

    void applyPrecon();

    void applyPressure();

    void advect();

    void rungeKutta3(double &x, double &y);


    void project();


    void scaleAdd(std::vector<double> &curr, std::vector<double> &a, std::vector<double> &b, double s);
};


static bool absCompare(double a, double b) {
    return (std::abs(a) < std::abs(b));
}


