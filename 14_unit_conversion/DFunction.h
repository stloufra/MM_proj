#ifndef DFUNCTION_H
#define DFUNCTION_H
#include <vector>
#include <cassert>

#pragma once

class DFunction
{
public:
    DFunction()=delete;
    DFunction(int r, int c)
    {
        rows_pr = r; //klasicke pocitani
        collums_pr = c;
        data.resize(rows_pr*collums_pr*velocities);
        collums = c;
        rows = r;
    }
    ~DFunction() {} 

    double &operator () (int i, int j, int v)
    {   
        assert(i >= 0 && i <= rows_pr && j >= 0 && j <= collums_pr && v <= velocities -1 && v >= 0);

        return data[i*collums_pr*velocities + j*velocities + v ];
    }

    const double &operator () (int i, int j, int v) const
    {   
         assert(i >= 0 && i <= rows_pr && j >= 0 && j <= collums_pr && v <= velocities -1 && v >= 0);

        return data[i*collums_pr*velocities + j*velocities + v ];
    }
    

    int rows;
    int collums;

private:
    std::vector<double> data;
    int rows_pr;
    int collums_pr;
    int const velocities=9;
};

#endif