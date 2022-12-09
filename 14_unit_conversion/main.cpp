#include <iostream>
#include <fstream>
#include "Mesher.h"
#include "Mesh.h"
#include "Solver.h"



    
int main()
{   
    //variables
    const int Nx = 200;
    const int Ny = 100;
    const double rho0=1;
    const double ux0=0.0;
    const double uy0=0;
    const double uw=0.15;
    const int iterations = 5000;
    const int plot_every=25;

    double Re = 400;
    

    Mesher mesh_rectangle(Ny,Nx);

    //object functions
    mesh_rectangle.obj_horizontal_line(mesh_rectangle.objects.collums); // -1, 0, 1, ..., Nx
    mesh_rectangle.obj_horizontal_line(-1);
    mesh_rectangle.obj_vertical_line(mesh_rectangle.objects.rows);
    mesh_rectangle.obj_vertical_line(-1);
    mesh_rectangle.obj_cylinder(Ny/2,Nx/4,Ny/5);
    mesh_rectangle.obj_rectangle(Ny/2,Nx*3/4,Ny/3,Ny/3);

    //meshing functions
    mesh_rectangle.meshing();
    mesh_rectangle.mesh_mowing_wall_horizontal(0, uw);

    //output mesh
    mesh_rectangle.output_VTK();

    //solver
    Solver solver(Ny,Nx,mesh_rectangle);
    solver.convert_to_lattice(0.2, 1000, 1, 1e-6);
    solver.initialization_eq(rho0, ux0, uy0);
    solver.output_VTK(0,plot_every);

    int k = 0;
    while(k<iterations) //err>=10e-4)
    {
        k++;
        solver.collision();
        solver.streaming();
        solver.bounce_back();
        solver.postpro();
        
        if(k%1000==0 && k!=0)
        {
            solver.Err();
            printf("err=%e ux_center=%e uy_center=%e rho_center=%e k=%d\n",solver.err,solver.ux(Ny/2,Nx/2),solver.uy(Ny/2,Nx/2),solver.rho(Ny/2,Nx/2), k);
        }

        if(k%plot_every==0)
        {
            solver.output_VTK(k,plot_every);
            std::cout<<k<<"\n";
        }
        
        
    }

    

    return 0;
}
/*
int main(){
Mesh Mesh(2,7);

    for(int j=-1; j < Mesh.rows+1; j++)
    {
        for(int i=-1; i < Mesh.collums+1; i++)
        {
            Mesh(j,i)=j+1;

        }

    }

    for(int j=-1; j < Mesh.rows+1; j++)
    {
        for(int i=-1; i < Mesh.collums+1; i++)
        {
            std::cout<<Mesh(j,i)<<" ";
        }
        std::cout<<"\n";
    }
return 0;
}
*/