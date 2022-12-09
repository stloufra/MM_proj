#ifndef MESHER_H
#define MESHER_H

#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include "Mesh.h"


class Mesher
{
public:
    Mesher()=delete;
    Mesher(int Ny_in, int Nx_in);

    //objects functions
    void obj_cylinder(double y_cylinder, double x_cylinder, double D_cylinder);
    void obj_rectangle(double y_rectangle, double x_rectangle, double lenght_rectangle, double height_rectangle);
    void obj_vertical_line(int Ny);
    void obj_horizontal_line(int Nx);

    //mesh functions
    void meshing();
    void mesh_mowing_wall_horizontal(int Ny, double uw_in);

    //help functions
    void output_VTK();
    double distance(double x1, double y1,double x2, double y2);
    
    ~Mesher();

    Mesh objects;
    Mesh mesh;

    double uw;
private:
    
};

#endif