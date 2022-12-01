#ifndef VARIABLES_H
#define VARIABLES_H
#include <vector>
#include <cassert>

#pragma once

class Variables
{
public:
    Variables()=delete;
    Variables(int r, int c)
    {
        rows_pr = r; //klasicke pocitani
        collums_pr = c;
        data.resize(rows_pr*collums_pr);
        collums = c;
        rows = r;
    }
    ~Variables() {} 

    double &operator () (int i, int j)
    {   
        assert(i >= 0 && i <= rows_pr && j >= 0 && j <= collums_pr );

        return data[i*collums_pr + j ];
    }

    const double &operator () (int i, int j) const
    {   
        assert(i >= 0 && i <= rows_pr && j >= 0 && j <= collums_pr );

        return data[i*collums_pr + j ];
    }
    

    int rows;
    int collums;

private:
    std::vector<double> data;
    int rows_pr;
    int collums_pr;
};


#endif