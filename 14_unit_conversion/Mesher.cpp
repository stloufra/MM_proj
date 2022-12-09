// OBJECTS - structured bolean values of solid
// 0 = solid | 1 = fluid
//
// MESH - structured int values of nodes types
//  0 0 0 0 0 0 0 0 0 0 0 
//    -------------------
//  0|9 4  4 4 4 4  4 8|0
//  0|5 1  1 1 1 1  1 3|0
//  0|5 1 13 2 2 12 1 3|0
//            ---
//  0|5 1  3|0 0|5 1  3|0
//            ---            
//  0|5 1 10 4 4 11 1 3|0
//  0|5 1 1  1 1  1 1 3|0
//  0|6 2 2  2 2  2 2 7|0
//    -------------------
//  0 0 0 0 0 0 0 0 0 0 0 



#include "Mesher.h"


Mesher::Mesher(int Ny_in, int Nx_in):
objects(Ny_in, Nx_in), //constructor
mesh(Ny_in, Nx_in) //constructor

{
  
  for(int j=-1; j < objects.rows+1; j++)
  {
    for(int i=-1; i < objects.collums+1; i++)
    {
      objects(j,i)=1;
    }

  }
 
  

  for(int j=-1; j < mesh.rows+1; j++)
  {
    for(int i=-1; i < mesh.collums+1; i++)
    {
      mesh(j,i)=1;
    }
  }

}


//--------------------------OBJECTS FUNCTIONS------------------------------------- 
void Mesher::obj_cylinder(double y_cylinder, double x_cylinder, double D_cylinder)
{
  for(int j=-1; j < objects.rows+1; j++)
  {
    for(int i=-1; i < objects.collums+1; i++)
    {
      if(distance(x_cylinder,y_cylinder,i,j) <= D_cylinder)
      {
        objects(j,i)=0;
      }
    }
  }
}

void Mesher::obj_vertical_line(int Ny)
{
  for(int i=-1; i < objects.collums+1; i++)
  {
    objects(Ny,i)=0;
  }
}

void Mesher::obj_horizontal_line(int Nx)
{
  for(int j=-1; j < objects.rows+1; j++)
  {
    objects(j,Nx)=0;
  }
}

void Mesher::obj_rectangle(double y_rectangle, double x_rectangle, double lenght_rectangle, double height_rectangle) 
{
  for(int j=-1; j < objects.rows+1; j++)
  {
    for(int i=-1; i < objects.collums+1; i++)
    {
      if( j < y_rectangle+height_rectangle/2 && j > y_rectangle-height_rectangle/2 && i > x_rectangle - lenght_rectangle/2 && i < x_rectangle + lenght_rectangle/2 )
      {
        objects(j,i)=0;
      }
    }
  }    
}

//--------------------------MESHING FUNCTIONS---------------------------------------


void Mesher::meshing()
{
  for(int j=-1; j < objects.rows+1; j++)
  {
    for(int i=-1; i < objects.collums+1; i++)
    {  
      //solid
      if (objects(j,i)==0)
      {
        mesh(j,i)=0;
      }
      
      //fluid
      if (objects(j,i)==1)
      {
        mesh(j,i)=1;
      }
      if(objects(j,i)==1)
      {
        
        //dolni
        if (objects(j-1,i) == 0) 
        {
          mesh(j,i)=2;
        }

        //prava
        if ( objects(j,i+1) == 0)
        {
          mesh(j,i)=3;
        }
        
        //horni
        if ( objects(j+1,i) == 0)
        {
          mesh(j,i)=4;
        }
        
        //leva
        if (objects(j,i-1) == 0)
        {
          mesh(j,i)=5;
        }

        //levy dolni vnnitrni
        if (objects(j-1,i-1) == 0 && objects(j,i-1) == 0 && objects(j-1,i) == 0)
        {
          mesh(j,i) = 6;
        }

        //pravy dolni vnnitrni
        if (objects(j-1,i+1) == 0 && objects(j,i+1) == 0 && objects(j-1,i) == 0)
        {
          mesh(j,i) = 7;
        }

        //pravy horni vnnitrni
        if (objects(j+1,i+1) == 0 && objects(j,i+1) == 0 && objects(j+1,i) == 0)
        {
          mesh(j,i) = 8;
        }

        //levy horni vnnitrni
        if (objects(j+1,i-1) == 0 && objects(j,i-1) == 0 && objects(j+1,i) == 0)
        {
          mesh(j,i) = 9;
        }

        //levy dolni vnejsi
        if (objects(j+1,i+1) == 0 && objects(j,i+1) == 1 && objects(j+1,i) == 1)
        {
          mesh(j,i) = 10;
        }

        //pravy dolni vnejsi
        if (objects(j+1,i-1) == 0 && objects(j,i-1) == 1 && objects(j+1,i) == 1)
        {
          mesh(j,i) = 11;
        }

        //pravy horni vnejsi
        if (objects(j-1,i-1) == 0 && objects(j,i-1) == 1 && objects(j-1,i) == 1)
        {
          mesh(j,i) = 12;
        }

        //levy horni vnejsi
        if (objects(j-1,i+1) == 0 && objects(j,i+1) == 1 && objects(j-1,i) == 1)
        {
          mesh(j,i) = 13;
        }
      }
    }
  }
  std::cout<<"Meshing undergone successfully."<<"\n";
}

void Mesher::mesh_mowing_wall_horizontal(int Ny, double uw_in)
{
  
  for(int i=-1; i < mesh.collums+1; i++)
  {  
  if(mesh(Ny,i)==2 || mesh(Ny,i)==6 || mesh(Ny,i)==7)
    {
      mesh(Ny,i)=22;
    }

    if(mesh(Ny,i)==4 || mesh(Ny,i)==9 || mesh(Ny,i)==8)
    {
      mesh(Ny,i)=24;
    }
  }

  uw=uw_in;

}

//--------------------------HELP FUNCTIONS---------------------------------------

double Mesher::distance(double x1, double y1,double x2, double y2)
{
    //simple distance calculation
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

void Mesher::output_VTK()
{

  int Nx2 = objects.collums +2;
  int Ny2 = objects.rows +2;

  std::ofstream out_file("results/Mesh.vtk");

  out_file << "# vtk DataFile Version 2.0\n" ;
  out_file << "LBE mesh\n" ;
  out_file << "ASCII\n";
  out_file << "DATASET STRUCTURED_POINTS\n" ;
  out_file << "DIMENSIONS "<<Nx2<<" "<<Ny2<<" 1\n" ;
  out_file << "ASPECT_RATIO 1 1 1\n";
  out_file << "ORIGIN 0 0 0\n";
  out_file << "POINT_DATA "<<Nx2*Ny2<<"\n";

  out_file << "SCALARS "<<"objects "<<"int 1\n";
  out_file << "LOOKUP_TABLE default\n";
  for(int j=-1; j < objects.rows+1; j++)
  {
    for(int i=-1; i < objects.collums+1; i++)
    { 
      out_file <<objects(j,i)<<"\n";
    }
  }

  out_file << "SCALARS "<<"mesh "<<"int 1\n";
  out_file << "LOOKUP_TABLE default\n";
  for(int j=-1; j < mesh.rows+1; j++)
  {
    for(int i=-1; i < mesh.collums+1; i++)
    {
      out_file <<mesh(j,i)<<"\n";
    }
  }
}

Mesher::~Mesher()
{

}
