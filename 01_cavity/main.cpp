#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// -----------------------------------------------------------------------------------------
// ------------------------ CONSTANTS AND DECLARATION --------------------------------------
// -----------------------------------------------------------------------------------------

const int Nx = 256 ; //even
const int Ny = 256 ; //even
const int Nx1 = Nx +1 ;
const int Ny1 = Ny +1 ;
const int L = Ny +1 ;
const int Nvel = 9 ; // number of possible descrete velocities
const int plot_every = 20;

const int cx_pos [Nvel] = { 0, 1, 0, -1, 0, 1, -1, -1,  1};
const int cy_pos [Nvel] = { 0, 0, 1, 0, -1, 1,  1, -1, -1};
const double weight [Nvel] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};

double rho0 = 1.0;
double ux0 = 0.0;
double uy0 = 0.0;
double uw = 0.1;
double Re = 400.0;
int k;

//declaration of arrays of variables

double rho [Nx][Ny]; 
double ux [Nx][Ny]; 
double uy [Nx][Ny]; 

//declaration of array of distribution function

double f [Nx] [Ny] [Nvel]; 
double f_post_col [Nx] [Ny] [Nvel]; 
double tau;
int Rvel[Nvel]={0, 3, 4, 1, 2, 7, 8, 5, 6 };
double ux0r [Nx1][Ny1];
double uy0r [Nx1][Ny1];

void initialization(void) ;
double f_equilibrium(double rho_eq, double ux_eq, double uy_eq, int k_eq);
void collision(void);
void streaming(void);
void bounce_back(void);
void postpro(void);
double Err(void);
double u0[Nx1][Ny1];
double v0[Nx1][Ny1];
void output(int k);


// -----------------------------------------------------------------------------------------
// ---------------------------------------- MAIN -------------------------------------------
// -----------------------------------------------------------------------------------------

int main()
{
    int X2, Y2;
    double err;
    X2=Nx/2; Y2=Ny/2;

    k=0;
    err=1.0;
    tau = 3*L*uw/Re +0.5;

    initialization();

    while(err>1.0e-6)
    {
        k++;
        collision();
        streaming();
        bounce_back();
        postpro();

        if(k%1000==0)
        {
            err=Err();
            printf("err=%e ux?center=%e uy_center=%e k=%d\n",err,ux[X2][Y2],uy[X2][Y2], k);
        }

        if(k%plot_every==0)
        {
            output(k);
        }

    }
    return 0;
}    

void initialization() 
{
    // initialization with equilibrium method
    int i, j ,k;

    for (i=0; i<Nx; i++)
    {
        for (j=0; j<Ny; j++)
        {
            rho [i][j] = rho0;
            ux [i][j] = ux0;
            uy [i][j] = uy0;

            for(k=0;k<Nvel;k++)
            {
                f[i][j][k] = f_equilibrium(rho [i][j], ux[i][j], uy[i][j], k);
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
}

void streaming()
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
}

void bounce_back()
{
    int i,j;

    //lower wall
    for(i=0;i<Nx; i++)
    {
        f[i][0][2]=f_post_col[i][0][4];
        f[i][0][5]=f_post_col[i][0][7];
        f[i][0][6]=f_post_col[i][0][8];
    }

    //left wall
    for(j=0;j<Ny; j++)
    {
        f[0][j][1]=f_post_col[0][j][3];
        f[0][j][5]=f_post_col[0][j][7];
        f[0][j][8]=f_post_col[0][j][6];
    }

    //right wall
    for(j=0;j<Ny; j++)
    {
        f[Nx][j][3]=f_post_col[Nx][j][1];
        f[Nx][j][7]=f_post_col[Nx][j][5];
        f[Nx][j][6]=f_post_col[Nx][j][8];
    }
   //upper wall
    for(i=0;i<Nx; i++)
    {
        f[i][Ny][4]=f_post_col[i][Ny][2];
        f[i][Ny][7]=f_post_col[i][Ny][5]+6*rho[i][Ny]*weight[7]*cx_pos[7]*uw;
        f[i][Ny][8]=f_post_col[i][Ny][6]+6*rho[i][Ny]*weight[8]*cx_pos[8]*uw;
    }

}
void postpro()
{
    int i, j;
    for (i=0; i<Nx; i++)
        {
            for (j=0; j<Ny; j++)
            { 
                rho[i][j] = f[i][j][0] + f[i][j][1] + f[i][j][2] + f[i][j][3] + f[i][j][4] + f[i][j][5] + f[i][j][6] + f[i][j][7] + f[i][j][8];
                ux[i][j] = (f[i][j][1] + f[i][j][5] + f[i][j][8] - f[i][j][3] - f[i][j][6] - f[i][j][7])/rho[i][j];
                ux[i][j] = (f[i][j][5] + f[i][j][6] + f[i][j][2] - f[i][j][7] - f[i][j][8] - f[i][j][4])/rho[i][j];
            }
        }
}



void output(int s)
{
    int i, j;

    std::string step =std::to_string(s/plot_every);

    std::ofstream out_file("LBE.csv."+step);

    out_file << "x coord, y coord, ux, uy, rho, z coord\n" ;
    
    for (i=0; i<Nx; i++)
        {
            for (j=0; j<Ny; j++)
            {
                out_file<<i<<", "<<j<<", "<<ux[i][j]<<", "<<uy[i][j]<<", "<<rho[i][j]<<", 1\n";
            }
        }

    out_file.close();


}

double Err(void) 
{
    int i, j;
    double er1, er2;
    er1=0; er2=0;

    for (i=1; i<Nx; i++)
    {
        for (j=0; j<Ny; j++)
        {
            er1=sqrt((ux[i][j]-ux0r[i][j])*(ux[i][j]-ux0r[i][j]) + (uy[i][j]-uy0r[i][j])*(ux[i][j]-uy0r[i][j]));
            er2=sqrt(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]);
            ux0r[i][j]=ux[i][j];
            uy0r[i][j]=uy[i][j];
            
        }
    }
    return er1/er2;
}