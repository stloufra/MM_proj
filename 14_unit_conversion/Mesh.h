#ifndef MESH_H
#define MESH_H
#include <vector>
#include <cassert>

#pragma once



class Mesh
// cislovani od -1 do N+1
{
public:

    Mesh() = delete;
    //{
    //}

    Mesh(int r, int c)
    {
        rows_pr = r + 2; //automaticke pridani okrajovych sten 1 ---- 1
        collums_pr = c + 2;
        data.resize(rows_pr*collums_pr);
        collums = c;
        rows = r;
    }
   

    int &operator () (int j, int i)
    {   
        assert(j >= -1 && j <= rows_pr-1 && i >= -1 && i < collums_pr-1 );

        return data[j*collums_pr + collums_pr + i+1];
    }

    const int &operator () (int i, int j) const
    {   
        assert(j >= -1 && j <= rows_pr-1 && i >= -1 && i < collums_pr-1 );

        return data[i*collums_pr* + collums_pr + j+1];
    }
    
    ~Mesh() {} 

    int rows;
    int collums;

private:
    std::vector<int> data;
    int rows_pr;
    int collums_pr;
    
};

#endif