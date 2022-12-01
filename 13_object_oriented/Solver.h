#ifndef SOLVER_H
#define SOLVER_H

#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include "Mesher.h"
#include "Variables.h"
#include "DFunction.h"

class Solver
{
public:
    Solver()=delete;
    Solver(int Nx_in, int Ny_in, Mesher mesh_in);
    void initialization_eq(double rho0, double ux0, double uy0, double Re);
    double f_equilibrium(double rho_eq, double ux_eq, double uy_eq, int k_eq);
    void collision();
    void streaming();
    void bounce_back();
    void postpro();
    void Err();
    void output_VTK(int s,int plot_every);
    ~Solver();

    double err = 1;

    Variables rho;
    Variables ux;
    Variables uy;

private:
    //Model parameters
    /*
    6    2    5
    \  |  /
    3 -- 0 -- 1 
    /  |  \
    7    4    8
    */

    const int Nvel = 9 ; // Q number of possible descrete velocities
    const int cx_pos[9] = { 0, 1, 0, -1, 0, 1, -1, -1,  1};
    const int cy_pos[9] = { 0, 0, 1, 0, -1, 1,  1, -1, -1};
    const double weight[9] ={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
    double tau;
    
    DFunction df;
    DFunction df_post;
    Mesher mesh;
    Variables ux0er;
    Variables uy0er;

    int Nx;
    int Ny;
    int s;

    double uw22=0.1;
    double uw24=0;

};

#endif