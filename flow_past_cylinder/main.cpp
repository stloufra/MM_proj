#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

// -----------------------------------------------------------------------------------------
// ------------------------ CONSTANTS AND DECLARATION --------------------------------------
// -----------------------------------------------------------------------------------------

const int Nx = 400 ;
const int Ny = 100 ;
const int Nt = 300 ;
const int plot_every = 25;

// cylinder
const double D_cylinder = 14 ;
const double x_cylinder = 100 ;
const double y_cylinder = 50 ;
bool objects [Nx][Ny];


//using D2Q9 model

const int Nvel = 9 ; // number of possible descrete velocities

const int cx_pos [9] = { 0, 1, 0, -1, 0, 1, -1, -1,  1};
const int cy_pos [9] = { 0, 0, 1, 0, -1, 1,  1, -1, -1};
const double weight [9] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};

//initial conditions

double rho0 = 1.0;
double ux0 = 1.0;
double uy0 = 0.0;
double const ny = 1.035e-6; //water

//declaration of arrays of variables

double rho [Nx][Ny]; 
double ux [Nx][Ny]; 
double uy [Nx][Ny]; 

//declaration of array of distribution function

double f [Nx] [Ny] [Nvel]; 
double f_post_col [Nx] [Ny] [Nvel]; 

// -----------------------------------------------------------------------------------------
// -------------------------------- HELPING FUNCTIONS --------------------------------------
// -----------------------------------------------------------------------------------------

double Re_f()
{
    //returns reynolds number
   ///double u = sqrt(ux0*ux0+uy0*uy0);
   //return D_cylinder*u/ny;
   return 0;
}

double tau_f()
{
    //returns relaxation time
    double u = sqrt(ux0*ux0+uy0*uy0);
    double Re = 300;//Re_f();
    return 3*D_cylinder*u/Re +0.5;
}

double distance(int x1, int y1,int x2, int y2)
{
    //simple distance calculation
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int boundary_objects()
{
    //structured bolean mesh values of existencion of objects
    int i, j, k;
     for (i=0; i<Nx; i++)
    {
        for (j=0; j<Ny; j++)
        {    
            if(distance(x_cylinder,y_cylinder,i,j)<D_cylinder)
            {
                objects[i][j]=true;
            }
            else
            {
                objects[i][j]=false;
            }
        }
     }
   
return 0;
}

// -----------------------------------------------------------------------------------------
// ----------------------------- EVOLUTION FUNCTIONS  --------------------------------------
// -----------------------------------------------------------------------------------------


double f_equilibrium(double rho_eq, double ux_eq, double uy_eq, int k_eq)
{
    double cu, u2;

    cu = cx_pos[k_eq]*ux_eq + cy_pos[k_eq]*uy_eq;
    u2 = ux_eq*ux_eq + uy_eq*uy_eq;

    return weight[k_eq]*rho_eq*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u2);                                    //equlibrium
}

int initialization_by_u() 
{
    // initialization with equilibrium method
    int i, j ,k;

    boundary_objects();

    for (i=0; i<Nx; i++)
    {
        for (j=0; j<Ny; j++)
        {
            
            ux [i][j] = ux0;
            uy [i][j] = uy0;

            rho [i][j] = rho0;
            for(k=0;k<Nvel;k++)
            {
                f[i][j][k] = f_equilibrium(rho [i][j], ux[i][j], uy[i][j], k);
                
                //boundary of object
                if(distance(x_cylinder,i,y_cylinder,j)<D_cylinder)
                {
                    f[i][j][k] = 0;
                }
            }

        }
    }    

return 0;
}

int collision(double tau) 
{
    int i, j ,k;
    double feq;
    

    for (i=0; i<Nx; i++)
    {
        for (j=0; j<Ny; j++)
        {    
             for(k=0;k<Nvel;k++)
            {
                feq = f_equilibrium(rho [i][j], ux[i][j], uy[i][j], k);
                f_post_col[i][j][k]= f[i][j][k] - (f[i][j][k]-feq)/tau;                             //collision
            }
        }
    }

    return 0;
}


int streaming()
{
    int i, j, id, jd, k;
    for (i=0; i<Nx; i++)
    {
        for (j=0; j<Ny; j++)
        {    
             for(k=0;k<Nvel;k++)
            {
                id = i - cx_pos[k];
                jd = j - cy_pos[k];

                if(id>=0 && id<=Nx && jd>=0 && jd<=Ny)
                {
                    f[i][j][k]=f_post_col[id][jd][k];                                               //streaming
                }

            }
        }
    }
    return 0;
}


int bounce_back()
{
    int i,j;

    //lower wall
    for(i=0;i<Nx; i++)
    {
        f[i][0][2]=f_post_col[i][0][4];
        f[i][0][5]=f_post_col[i][0][8];
        f[i][0][6]=f_post_col[i][0][7];
    }

    //upper wall
    for(i=0;i<Nx; i++)
    {
        f[i][Ny][4]=f_post_col[i][Ny][2];
        f[i][Ny][8]=f_post_col[i][Ny][5];
        f[i][Ny][7]=f_post_col[i][Ny][6];
    }

    //boundary objects
    for (i=0; i<Nx; i++)
    {
        for (j=0; j<Ny; j++)
        { 
            if(objects[i][j])
            {
                f[i][j][3]=f_post_col[i][j][1];
                f[i][j][4]=f_post_col[i][j][2];
                f[i][j][1]=f_post_col[i][j][3];
                f[i][j][2]=f_post_col[i][j][4];
                f[i][j][7]=f_post_col[i][j][5];
                f[i][j][8]=f_post_col[i][j][6];
                f[i][j][5]=f_post_col[i][j][7];
                f[i][j][6]=f_post_col[i][j][8];

                
            }          
        }
    }
    //Zou-He boundary condition 

    /*for (j=0, j<Ny, j++)
    {
        f[Ny][j][6]
    }*/
    return 0;
}

// -----------------------------------------------------------------------------------------
// ----------------------------- POST PROCESSING FUNCTIONS ---------------------------------
// -----------------------------------------------------------------------------------------

int postpro()
{
    int i, j;
    for (i=0; i<Nx; i++)
        {
            for (j=0; j<Ny; j++)
            { 
                rho[i][j] = f[i][j][0] + f[i][j][1] + f[i][j][2] + f[i][j][3] + f[i][j][4] + f[i][j][5] + f[i][j][6] + f[i][j][7] + f[i][j][8];
                ux[i][j] = (f[i][j][1] + f[i][j][5] + f[i][j][8] - f[i][j][3] - f[i][j][6] - f[i][j][7])/rho[i][j];
                ux[i][j] = (f[i][j][5] + f[i][j][6] + f[i][j][2] - f[i][j][7] - f[i][j][8] - f[i][j][4])/rho[i][j];

                //inside cylinder 

                if(objects[i][j])
                {
                    ux[i][j] = 0.0;
                    uy[i][j] = 0.0;
                    rho[i][j]= 0.0;
                }
            }
        }
return 0;
}

void test(int t)
{
    std::cout << "\n f > "<< f[200][20][1] << "\n";
    std::cout << "\n rho > "<< rho[200][20] << "\n";
    std::cout << "\n ux > "<< ux[400][100] << "\n";

    std::cout << "Time > "<< t << "\n";
}

void output(int s)
{
    int i, j;

    std::string step =std::to_string(s/plot_every);

    std::ofstream out_file("LBE.csv."+step);

    out_file << "x coord, y coord, ux, uy, rho, obj\n" ;
    
    for (i=0; i<Nx; i++)
        {
            for (j=0; j<Ny; j++)
            {
                out_file<<i<<", "<<j<<" ,"<<ux[i][j]<<" ,"<<uy[i][j]<<" ,"<<rho[i][j]<<" ,"<<objects[i][j]<<"\n";
            }
        }

    out_file.close();


}

// -----------------------------------------------------------------------------------------
// ---------------------------------------- MAIN -------------------------------------------
// -----------------------------------------------------------------------------------------

int main()
{
    int t,s=0;
    double tau = tau_f();
    initialization_by_u();

    for (t=0;t<Nt;t++)
    {
        collision(tau);
        streaming();
        bounce_back();
        postpro();
        test(t);
        if(t==s|t==Nt)
        {
            output(s);
            s+=plot_every;
        }
    }


    
    return 0;
}
