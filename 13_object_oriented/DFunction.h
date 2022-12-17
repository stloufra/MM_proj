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

    double &operator () (int y, int x, int v)
    {   
        assert(y >= 0 && y < rows_pr && x >= 0 && x < collums_pr && v < velocities  && v >= 0);

        return data[y*collums_pr*velocities + x*velocities + v ];
    }

    const double &operator () (int y, int x, int v) const
    {   
        assert(y >= 0 && y < rows_pr && x >= 0 && x < collums_pr && v < velocities  && v >= 0);

        return data[y*collums_pr*velocities + x*velocities + v ];
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