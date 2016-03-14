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

    /*
     * Update simulation one timestep
     */
    void update();

    /*
     * Get timestep
     */
    double getTimestep();

    /*
     * Get pointer to image
     */
    unsigned char *getImagePtr();

    /*
     * Update the image with the current simulation data
     */
    void updateImage();


private:


    /*
     * Constants
     */

    const double _GRAVITY = 9.82; // Gravity
    const double _RHO = 0.1; // Density of fluid
    const double _DENSITYSOOT = 0.25;
    const double _DENSITYAIR = 0.1;
    const double TAMB = 273; // Ambient temperature
    const double KDISS = 0.001; // Smoke dissipation factor
    const int _ITERLIMIT = 600; // Iteration limit in the PCG solver
//    const double PI = 3.14159265;
    const double _TOL = 10e-5; // Tolerance in the PCG solver

    // MIC constants
    double _tau_mic = 0.97;
    double _sigma_mic = 0.25;


    /*
     * Member variables
     */

//    int pulse_time;

    unsigned char *_image; // The image array

    // Number of cells in the x and y direction. Right now these have to be the same.
    unsigned int _nx;
    unsigned int _ny;

    unsigned int length; // Length of the vectors containing the cells
    double _dxy; // Distance between cells
    double _dt; // Timestep

    // Vectors storing the coefficient matrix
    std::vector<double> _Adiag;
    std::vector<double> _Aplusi;
    std::vector<double> _Aplusj;

    std::vector<double> _rhs; // Negative divergence vector, (right hand side)
    std::vector<double> _p; // Pressure vector
    std::vector<double> _s; // Temporary vector for the pcg algorithm
    std::vector<double> _z; // Temporary vector for the pcg algorithm
    std::vector<double> _precon; // Vector storing the preconditioner
    // Vectors storing the velocities in the x and y directions
    std::vector<double> _u;
    std::vector<double> _v;
    // and temporary variables for them
    std::vector<double> _un;
    std::vector<double> _vn;

    // Vector storing the smoke density/concentration in every cell
    std::vector<double> _d;
    // and its temporary variable
    std::vector<double> _dn;
    // Vector storing the temperature in every cell
    std::vector<double> _T;
    // Temporary variable
    std::vector<double> _Tn;


    /*
     * Functions
     */

    inline unsigned int getIdx(int i, int j, int width) {
        return i + j * width;
    }


    /*
     * Perform linear interpolation
     */
    double lerp(double a, double b, double x);


    /*
     * Perform bilinear interpolation
     */
    double lerp2(double x, double y, double ox, double oy, int w, int h, const std::vector<double> &quantity);


    //double lerp3(int x, int y, int z, float ox, float oy, float oz, int w, int h,
    //             std::unique_ptr<std::vector<double>> quantity);

    /*
     * Perform cubic interpolation
     */
    double cerp(double a, double b, double c, double d, double x);

    /*
     * Perform bicubic interpolation
     */
    double cerp2(double x, double y, int w, int h, double ox, double oy, std::vector<double> &quantity);


    /*
     * Wrapper for adding a source of smoke in the simulation
     */
    void addInFlow();

    /*
     * Add a source of fluids in the simulation
     */
    void addInFlow(double x0, double y0, double x1, double y1, int w, int h, double ox, double oy, double dxy,
                   double value,
                   std::vector<double> &src);


    /*
     * Smoke specific buoyancy calculations
     */
    void applyBuoyancy();

    /*
     * Calculate the negative divergence in each fluid cell
     */
    void buildRhs();

    /*
     * Build the coefficient matrix for the pressure equations
     */
    void buildPressureMatrix();

    /*
     * Build the MIC preconditioner
     */
    void buildPrecon();

    /*
     * Function to apply the preconditioner to the current pressure vector
     */
    void applyPrecon();

    /*
     * Update the velocities with the calculated pressure gradient
     */
    void applyPressure();

    /*
     * Self advect the velocities and advect the temperature and smoke
     */
    void advect();

    /*
     * Third order Runge-Kutta
     */
    void rungeKutta3(double &x, double &y);

    /*
     * Projection step
     */
    void project();


    /*
     * Utility function to add a scaled vector to another vector
     */
    void scaleAdd(std::vector<double> &curr, std::vector<double> &a, std::vector<double> &b, double s);
};





