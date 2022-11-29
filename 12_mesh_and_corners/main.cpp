
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

// -----------------------------------------------------------------------------------------
// ------------------------ CONSTANTS AND DECLARATION --------------------------------------
// -----------------------------------------------------------------------------------------

// Lattice and initial parameters 

const int Nx = 400 ; //even
const int Ny = 100 ; //even
const double rho0 = 1.0;
const double ux0 = 0.0;
const double uy0 = 0.0;
const double uw24 = 0.1;
const double uw22 = -0.1;
const double Re = 400;
const double ny = 1;

//Cylinder
bool cylinder_existence =0;
const double D_cylinder = 12 ;
const double x_cylinder1 = Nx/4 ;
const double y_cylinder1 = Ny/2 ;

//help

const int plot_every = 250; //25 for 5000

//mesh

const int Nx2 = Nx +2 ;
const int Ny2 = Ny +2 ;
bool objects [Ny2][Nx2];
int mesh [Ny2][Nx2];

//Model parameters

/*

6    2    5
  \  |  /
3 -- 0 -- 1 
  /  |  \
7    4    8

*/
const int Nvel = 9 ; // Q number of possible descrete velocities
const int cx_pos[Nvel] = { 0, 1, 0, -1, 0, 1, -1, -1,  1};
const int cy_pos[Nvel] = { 0, 0, 1, 0, -1, 1,  1, -1, -1};
const double weight[Nvel] ={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double tau;

//declaration of array of distribution function

double f[Ny2][Nx2][Nvel]; 
double f_post_col[Ny2][Nx2][Nvel]; 

//declaration of arrays of variables

double rho[Ny2][Nx2]; 
double ux[Ny2][Nx2]; 
double uy[Ny2][Nx2]; 

double ux0er[Ny2][Nx2];
double uy0er[Ny2][Nx2];

//functions 

void initialization_equilibrium(void) ;
double f_equilibrium(double rho_eq, double ux_eq, double uy_eq, int k_eq);
void collision(void);
void streaming(void);
void postpro(void);
void bounce_back(void);
double Err(void);

//help functions

double distance(int x1, int y1,int x2, int y2);
void boundary_objects(void);
void mesh_init(void);
void mesh_moving_walls_horizontal(int Ny_h);
void output_VTK_Variables(int s);
void output_VTK_StructPoint();
void inlet_outlet_BC(void);


// -----------------------------------------------------------------------------------------
// ---------------------------------------- MAIN -------------------------------------------
// -----------------------------------------------------------------------------------------

int main()
{
    int k,X2, Y2;
    double err;
    X2=Nx2/2; Y2=Ny2/2;
    k=0;
    err=1.0;
    tau = 3*Nx*uw24/Re +0.5;
    //tau = ny/3+0.5;
    
    //mesh operations
    boundary_objects();
    mesh_init();
    mesh_moving_walls_horizontal(1);
    mesh_moving_walls_horizontal(Ny);

    output_VTK_StructPoint();


    std::cout<<"Meshing succesfull.\n";
    
    //initialization
    initialization_equilibrium();

    output_VTK_Variables(k); 
    
    std::cout<<"Initialization succesfull.\n";

    while(k<5000) //err>=10e-4)
    {
        k++;
        collision();
        streaming();
        bounce_back();
        //inlet_outlet_BC();
        postpro();
        err=Err();
        if(k%1000==0)
        {
            printf("err=%e ux_center=%e uy_center=%e rho_center=%e k=%d\n",err,ux[Y2][X2],uy[Y2][X2],rho[Y2][X2], k);
        }

        if(k%plot_every==0)
        {
            output_VTK_Variables(k);
            std::cout<<k<<"\n";
        }
        
    }   
    output_VTK_Variables(k);
    
    return 0;
}    


// -----------------------------------------------------------------------------------------
// ----------------------------- EVOLUTION FUNCTIONS  --------------------------------------
// -----------------------------------------------------------------------------------------

void initialization_equilibrium() 
{
    // initialization with equilibrium method
    int j, i, k;

    for (j=0; j<Ny2; j++)
    {
        for (i=0; i<Nx2; i++)
        {
            if(objects[j][i]==1)
            {
                rho [j][i] = rho0;
                ux [j][i] = ux0;
                uy [j][i] = uy0;

                for(k=0;k<Nvel;k++)
                {
                    f[j][i][k] = f_equilibrium(rho [j][i], ux[j][i], uy[j][i], k);
                    
                }
            }
            if(objects[j][i]==0)
            {
                rho [j][i] = 0;
                ux [j][i] = 0;
                uy [j][i] = 0;
                for(k=0;k<Nvel;k++)
                {
                    f[j][i][k] =0;
                }
            }
        }
    }    
}

double f_equilibrium(double rho_eq, double ux_eq, double uy_eq, int k_eq)
{
    double cu, u2;

    cu = cx_pos[k_eq]*ux_eq + cy_pos[k_eq]*uy_eq;
    u2 = ux_eq*ux_eq + uy_eq*uy_eq;

    return weight[k_eq]*rho_eq*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u2);                                    //equlibrium
}

void collision() 
{
    int i, j ,k;
    double feq;
    

    for (j=0; j<Ny2; j++)
    {
        for (i=0; i<Nx2; i++)
        {   
            if(objects[j][i]==1)
            {
                for(k=0;k<Nvel;k++)
                {
                    feq = f_equilibrium(rho [j][i], ux[j][i], uy[j][i], k);
                    f_post_col[j][i][k]= f[j][i][k] - (f[j][i][k]-feq)/tau;                             //collision
                }
            }
        }
    }
}

void streaming()
{
    int j, i, jd, id, k;
    for (j=0; j<Ny2 ; j++)
    {
        for (i=0; i<Nx2; i++)
        {    
            for(k=0;k<Nvel;k++)
            {
                if(objects[j][i]==1)
                {
                    jd = j - cy_pos[k];
                    id = i - cx_pos[k];
                    
                    f[j][i][k]=f_post_col[jd][id][k];                                               //streaming
                    
                }

                if(objects[j][i]==0)
                {
                   f[j][i][k]=0;
                } 
            }
            
        }
    }
}

void bounce_back()
{
    
    for (int j=0; j<Ny2; j++)
    {
        for (int i=0; i<Nx2; i++)
        {
            //rovne stacionarni steny
            if(mesh[j][i]==2)
            {
                f[j][i][2]=f_post_col[j][i][4];
                f[j][i][5]=f_post_col[j][i][7];
                f[j][i][6]=f_post_col[j][i][8];
            }

            if(mesh[j][i]==3)
            {
                f[j][i][3]=f_post_col[j][i][1];
                f[j][i][6]=f_post_col[j][i][8];
                f[j][i][7]=f_post_col[j][i][5];
            }

            if(mesh[j][i]==4)
            {
                f[j][i][4]=f_post_col[j][i][2];
                f[j][i][7]=f_post_col[j][i][5];
                f[j][i][8]=f_post_col[j][i][6];
            }

            if(mesh[j][i]==5)
            {
                f[j][i][1]=f_post_col[j][i][3];
                f[j][i][5]=f_post_col[j][i][7];
                f[j][i][8]=f_post_col[j][i][6];
            }

            //rohy vnitrni
            if(mesh[j][i]==6)
            {
                f[j][i][2]=f_post_col[j][i][4];
                f[j][i][1]=f_post_col[j][i][3];
                f[j][i][5]=f_post_col[j][i][7];
            }

            if(mesh[j][i]==7)
            {
                f[j][i][2]=f_post_col[j][i][4];
                f[j][i][3]=f_post_col[j][i][1];
                f[j][i][6]=f_post_col[j][i][8];
            }

            if(mesh[j][i]==8)
            {
                f[j][i][4]=f_post_col[j][i][2];
                f[j][i][3]=f_post_col[j][i][1];
                f[j][i][7]=f_post_col[j][i][5];
            }

            if(mesh[j][i]==9)
            {
                f[j][i][4]=f_post_col[j][i][2];
                f[j][i][1]=f_post_col[j][i][3];
                f[j][i][8]=f_post_col[j][i][6];
            }

            //rohy vnejsi
            if(mesh[j][i]==10)
            {
                f[j][i][5]=f_post_col[j][i][7];
            }

            if(mesh[j][i]==11)
            {
                f[j][i][8]=f_post_col[j][i][6];
            }

            if(mesh[j][i]==12)
            {
                f[j][i][5]=f_post_col[j][i][7];
                
            }

            if(mesh[j][i]==13)
            {
                f[j][i][6]=f_post_col[j][i][8];
            }

            //rovne pohyblive steny
            if(mesh[j][i]==22)
            {
                f[j][i][2]=f_post_col[j][i][4];
                f[j][i][5]=f_post_col[j][i][7]+6*rho[j][i]*weight[5]*cx_pos[5]*uw22;
                f[j][i][6]=f_post_col[j][i][8]+6*rho[j][i]*weight[6]*cx_pos[6]*uw22; 
            }

            if(mesh[j][i]==24)
            {
                f[j][i][4]=f_post_col[j][i][2];
                f[j][i][7]=f_post_col[j][i][5]+6*rho[j][i]*weight[7]*cx_pos[7]*uw24;
                f[j][i][8]=f_post_col[j][i][6]+6*rho[j][i]*weight[8]*cx_pos[8]*uw24;
            }
        }
    }

   
    
}
void postpro()
{
    int i, j;
    for (j=0; j<Ny2; j++)
    {
        for (i=0; i<Nx2; i++)
        { 
            if(objects[j][i]==1)
            {
                rho[j][i]=f[j][i][0]+f[j][i][1]+f[j][i][2]+f[j][i][3]+f[j][i][4]+f[j][i][5]+f[j][i][6]+f[j][i][7]+f[j][i][8];
                ux[j][i]=(f[j][i][1]+f[j][i][5]+f[j][i][8]-f[j][i][3]-f[j][i][6]-f[j][i][7])/rho[j][i];
                uy[j][i]=(f[j][i][5]+f[j][i][6]+f[j][i][2]-f[j][i][7]-f[j][i][8]-f[j][i][4])/rho[j][i];      
            }
        }

    }

   

}

double Err() 
{
    int j, i;
    double er1, er2;
    er1=0.0; er2=0.0;

    for (j=1; j<=Ny; j++)
    {
        for (i=0; i<=Nx; i++)
        {
            er1+=sqrt((ux[j][i]-ux0er[j][i])*(ux[j][i]-ux0er[j][i])+(uy[j][i]-uy0er[j][i])*(uy[j][i]-uy0er[j][i]));
            er2+=sqrt(ux[j][i]*ux[j][i]+uy[j][i]*uy[j][i]);
            ux0er[j][i]=ux[j][i];
            uy0er[j][i]=uy[j][i];

           
        }
    }
    return er1/er2;
}



// -----------------------------------------------------------------------------------------
// -------------------------------- HELPING FUNCTIONS --------------------------------------
// -----------------------------------------------------------------------------------------

void boundary_objects()
{
    //structured bolean mesh values of solid
    // 0 = solid | 1 = fluid
    int i, j, k;
    for (j=1; j<=Ny; j++)
    {
        for (i=1; i<=Nx; i++)
        {   
            
                objects[j][i]=1;
        }
    }

    if(cylinder_existence)
    {
        for (j=1; j<=Ny; j++)
        {
            for (i=1; i<=Nx; i++)
            {   
                if(distance(x_cylinder1,y_cylinder1,i,j) <= D_cylinder)
                {
                    objects[j][i]=0;
                }
            }
        }
    }

    for (i = 0; i < Nx2; i++ )
    {
        objects[0][i] = 0;
        objects[Ny2][i] = 0;
    }

    for (i = 0; i < Ny2; i++ )
    {
        objects[i][0] = 0;
        objects[i][Nx2] = 0;
    }
}

void mesh_init()
{
    /*

    0 0 0 0 0 0 0 0 0 0 0 
      -------------------
    0|9 4  4 4 4 4  4 8|0
    0|5 1  1 1 1 1  1 3|0
    0|5 1 13 2 2 12 1 3|0
             ---
    0|5 1  3|0 0|5 1  3|0
             ---            
    0|5 1 10 4 4 11 1 3|0
    0|5 1 1  1 1  1 1 3|0
    0|6 2 2  2 2  2 2 7|0
      -------------------
    0 0 0 0 0 0 0 0 0 0 0 
        
    */
   for(int j = 0 ; j < Ny2; j++)
   {    
        for(int i = 0; i < Nx2; i++)
        {   
            //solid
            if (objects[j][i]==0)
            {
                mesh[j][i]=0;
            }
            
            //fluid
            if (objects[j][i]==1)
            {
                mesh[j][i]=1;
            }
            if(objects[j][i]==1)
            {
            
                //dolni
                if (objects[j-1][i] == 0) 
                {
                    mesh[j][i]=2;
                }

                //prava
                if ( objects[j][i+1] == 0)
                {
                    mesh[j][i]=3;
                }
                
                //horni
                if ( objects[j+1][i] == 0)
                {
                    mesh[j][i]=4;
                }
                
                //leva
                if (objects[j][i-1] == 0)
                {
                    mesh[j][i]=5;
                }

                //levy dolni vnnitrni
                if (objects[j-1][i-1] == 0 && objects[j][i-1] == 0 && objects[j-1][i] == 0)
                {
                    mesh[j][i] = 6;
                }

                //pravy dolni vnnitrni
                if (objects[j-1][i+1] == 0 && objects[j][i+1] == 0 && objects[j-1][i] == 0)
                {
                    mesh[j][i] = 7;
                }

                //pravy horni vnnitrni
                if (objects[j+1][i+1] == 0 && objects[j][i+1] == 0 && objects[j+1][i] == 0)
                {
                    mesh[j][i] = 8;
                }

                //levy horni vnnitrni
                if (objects[j+1][i-1] == 0 && objects[j][i-1] == 0 && objects[j+1][i] == 0)
                {
                    mesh[j][i] = 9;
                }

                //levy dolni vnejsi
                if (objects[j+1][i+1] == 0 && objects[j][i+1] == 1 && objects[j+1][i] == 1)
                {
                    mesh[j][i] = 10;
                }

                //pravy dolni vnejsi
                if (objects[j+1][i-1] == 0 && objects[j][i-1] == 1 && objects[j+1][i] == 1)
                {
                    mesh[j][i] = 11;
                }

                //pravy horni vnejsi
                if (objects[j-1][i-1] == 0 && objects[j][i-1] == 1 && objects[j-1][i] == 1)
                {
                    mesh[j][i] = 12;
                }

                //levy horni vnejsi
                if (objects[j-1][i+1] == 0 && objects[j][i+1] == 1 && objects[j-1][i] == 1)
                {
                    mesh[j][i] = 13;
                }

                
            }

        }
   }
}

/*
void inlet_outlet_BC(void)
{
    int j;
    for (j=0; j<Ny2; j++)
    {
        f[j][0][1]=f[j][Nx][1];
        f[j][0][5]=f[j][Nx][5];
        f[j][0][8]=f[j][Nx][8];

        f[j][Nx][6]=f[j][0][6];
        f[j][Nx][7]=f[j][0][7];
        f[j][Nx][3]=f[j][0][3];
    }
}
*/
void mesh_moving_walls_horizontal(int Ny_h)
{
    for(int i = 1; i<=Nx; i++)
    {
        if(mesh[Ny_h][i]==2)
        {
            mesh[Ny_h][i]=22;
        }

        if(mesh[Ny_h][i]==4)
        {
            mesh[Ny_h][i]=24;
        }
    }
}

double distance(int x1, int y1,int x2, int y2)
{
    //simple distance calculation
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}


void output_VTK_Variables(int s)
{
       
    std::string step =std::to_string(s/plot_every);

    std::ofstream out_file("results/LBE."+step+".vtk");

    out_file << "# vtk DataFile Version 2.0\n" ;
    out_file << "LBE mesh\n" ;
    out_file << "ASCII\n";
    out_file << "DATASET STRUCTURED_POINTS\n" ;
    out_file << "DIMENSIONS "<<Nx2<<" "<<Ny2<<" 1\n" ;
    out_file << "ASPECT_RATIO 1 1 1\n";
    out_file << "ORIGIN 0 0 0\n";
    out_file << "POINT_DATA "<<Nx2*Ny2<<"\n";

    out_file << "SCALARS "<<"rho "<<"double 1\n";
    out_file << "LOOKUP_TABLE default\n";
    for (int j=0; j<Ny2; j++)
    {
        for (int i=0; i<Nx2; i++)
        {
            out_file <<rho[j][i]<<"\n";
        }
    }
    out_file << "SCALARS "<<"objects "<<"double 1\n";
    out_file << "LOOKUP_TABLE default\n";
    for (int j=0; j<Ny2; j++)
    {
        for (int i=0; i<Nx2; i++)
        {
            out_file <<objects[j][i]<<"\n";
        }
    }
   
    out_file << "VECTORS "<<"U "<<"double\n";
    for (int j=0; j<Ny2; j++)
    {
        for (int i=0; i<Nx2; i++)
        {
            out_file <<ux[j][i]<<" "<<uy[j][i]<<" 0\n";
        }
    }


    out_file.close();


}

void output_VTK_StructPoint()
{
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
    for (int j=0; j<Ny2; j++)
    {
        for (int i=0; i<Nx2; i++)
        {
            out_file <<objects[j][i]<<"\n";
        }
    }

    out_file << "SCALARS "<<"mesh "<<"int 1\n";
    out_file << "LOOKUP_TABLE default\n";
    for (int j=0; j<Ny2; j++)
    {
        for (int i=0; i<Nx2; i++)
        {
            out_file <<mesh[j][i]<<"\n";
        }
    }
}
