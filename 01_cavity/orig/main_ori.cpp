//=========================================================
//---------------------------------------------------------
// ----- Header file of the D2Q9 model -----
//---------------------------------------------------------
//File name: D2Q9.h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#define Nx 256 // number of cells in the x-direction
#define Ny 256 // number of cells in the y-direction
#define Nx1 (Nx+1)
#define Ny1 (Ny+1)
#define L (Ny+1) // width of the cavity
#define Q 9 // number of discrete velocities
#define rho0 1.0 // initial density
#define ux0 0.0 // initial velocity component in x direction
#define uy0 0.0 // initial velocity component in y direction
#define uw 0.1
#define Re 400.0

int cx[Q]={0, 1, 0, -1, 0, 1, -1, -1, 1};
int cy[Q]={0, 0, 1, 0, -1, 1, 1, -1, -1};

double f[Ny1][Nx1][Q]; //array of the distribution functions (DFs)
double f_post[Ny1][Nx1][Q]; // array of the post-collision DFs

double rho[Ny1][Nx1];
double ux[Ny1][Nx1];
double uy[Ny1][Nx1];// arrays of fluid density and velocity

double u0[Ny1][Nx1];
double v0[Ny1][Nx1];

double tau; // relaxation time for BGK model

double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36}; // the weights in the EDF

//int rc[Q]={0,3,4,1,2,7,8,5,6}; // index of reversed velocity

void Init_Eq(void); //Initialization
double feq(double RHO, double U, double V, int k);// Equilibrium distribution function
void Coll_BGK(void); // BGK collision
void Streaming(void); // Streaming
void Den_Vel(void); // Fluid variables
void Bounce_back(void); // Bounce-back boundary condition
double Err(void); // Difference in velocity field
void Data_Output(void); // Output simulation data

//=========================================================
//=========================================================
int main()
{
    int k,M2,N2;
    double err;
    M2=Ny/2; N2=Nx/2;
    k=0;
    err=1.0;
    tau=3*L*uw/Re+0.5; // relaxation time for BGK
    
    Init_Eq();
    while(err>1.0e-3)
    {
        k++;
        Coll_BGK(); //BGK collision
        Streaming(); // Streaming
        Bounce_back(); // No-Slip boundary condition
        Den_Vel(); // Fluid variables
        if(k%1000==0)
        {
            err=Err(); // Velocity differences between two successive 1000 steps
            printf("err=%e ux_center=%e uy_center=%e k=%d\n",err,ux[M2][N2],uy[M2][N2], k); // Display some results
        }
    }
    Data_Output(); // Output simulation data
    return 0;
}
//=========================================================
//-------------------------------------------------------------------
// Subroutine: initialization with the equilibrium method
//------------------------------------------------------------------
//
void Init_Eq()
{
    int j, i, k;
    for (j=0;j<=Ny;j++) for(i=0;i<=Nx;i++)
    {
        rho[j][i]=rho0;
        ux[j][i]=ux0;
        uy[j][i]=uy0;
        for(k=0;k<Q;k++)
        f[j][i][k]=feq(rho[j][i],ux[j][i],uy[j][i],k);
    }
}
//========================================================
//=========================================================
//-----------------------------------------------------------------
// Subroutine: calculation the equilibrium distribution
//----------------------------------------------------------------
//
double feq(double RHO, double U, double V, int k)
{
    double cu, U2;
    cu=cx[k]*U+cy[k]*V; // c k*u
    U2=U*U+V*V; // u*u;
    return w[k]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}
//=========================================================
//=========================================================
//---------------------------------------------------------
// Subroutine: BGK collision
//---------------------------------------------------------
void Coll_BGK()
{
    int j, i, k;
    double FEQ;
    for (j=0;j<=Ny;j++) for(i=0;i<=Nx;i++) for(k=0;k<Q;k++)
    {
        FEQ=feq(rho[j][i],ux[j][i],uy[j][i],k); // EDF
        f_post[j][i][k] = f[j][i][k]-(f[j][i][k]-FEQ)/tau;
        // Post-collision DFs
    }
}
//=========================================================
//=========================================================
//---------------------------------------------------------
// Subroutine: Streaming
//---------------------------------------------------------
void Streaming()
{
    int j, i, jd, id, k;
    for (j=0;j<=Ny;j++) for(i=0;i<=Nx;i++) for(k=0;k<Q;k++)
    {
        jd=j-cy[k]; id=i-cx[k]; // upwind node
        if(jd>=0 && jd<=Ny && id>=0 && id<=Nx) // fluid node
        f[j][i][k]=f_post[jd][id][k]; // streaming
    }
}
//=========================================================
//=========================================================
//---------------------------------------------------------
// Subroutine: Bounce-back scheme
//---------------------------------------------------------
void Bounce_back()
{
    int i,j;
    // j=Ny: top plate
     for(i=0;i<=Nx; i++)
    {
        f[Ny][i][4]=f[Ny][i][2];
        f[Ny][i][7]=f[Ny][i][5]+6*rho[Ny][i]*w[7]*cx[7]*uw;
        f[Ny][i][8]=f[Ny][i][6]+6*rho[Ny][i]*w[8]*cx[8]*uw;
    }
    // j=0: bottom plate
     for(i=0;i<=Nx;i++)
    {
        f[0][i][2]=f_post[0][i][4];
        f[0][i][5]=f_post[0][i][7];
        f[0][i][6]=f_post[0][i][8];
    }
    // i=0: left wall
    for(j=0;j<=Ny;j++)
    {
        f[j][0][1]=f_post[j][0][3];
        f[j][0][5]=f_post[j][0][7];
        f[j][0][8]=f_post[j][0][6];
    }
    // i=Nx: right wall
    for(j=0;j<=Ny;j++)
    {
        f[j][Nx][3]=f_post[j][Nx][1];
        f[j][Nx][7]=f_post[j][Nx][5];
        f[j][Nx][6]=f_post[j][Nx][8];
    }
}
//=========================================================
//=========================================================
//------------------------------------------------------------
// Subroutine: Fluid variables (density and velocity)
//------------------------------------------------------------
void Den_Vel()
{
    int j, i;
    for(j=0;j<=Ny;j++) for(i=0;i<=Nx;i++)
    {
        rho[j][i]=f[j][i][0]+f[j][i][1]+f[j][i][2]+f[j][i][3]+f[j][i][4]+f[j][i][5]+f[j][i][6]+f[j][i][7]+f[j][i][8];
        ux[j][i]=(f[j][i][1]+f[j][i][5]+f[j][i][8]-f[j][i][3]-f[j][i][6]-f[j][i][7])/rho[j][i];
        uy[j][i]=(f[j][i][5]+f[j][i][6]+f[j][i][2]-f[j][i][7]-f[j][i][8]-f[j][i][4])/rho[j][i];
    }
}
//=========================================================
double Err() // Calculating the relative difference in velocity between two steps
{
    int j, i;
    double e1,e2;
    e1=e2=0.0;
    for(j=1;j<Ny;j++) for(i=0;i<Nx;i++)
    {
        e1+=sqrt((ux[j][i]-u0[j][i])*(ux[j][i]-u0[j][i])+(uy[j][i]-v0[j][i])*(uy[j][i]-v0[j][i]));
        e2+=sqrt(ux[j][i]*ux[j][i]+uy[j][i]*uy[j][i]);
        u0[j][i]=ux[j][i];v0[j][i]=uy[j][i];
    }
    return e1/e2;
}
void Data_Output() // Output data
{
    int i, j;

    std::ofstream out_file("LBE.csv");

    out_file << "x coord, y coord, ux, uy, rho, z coord\n" ;
    
    for (i=0; i<Nx; i++)
        {
            for (j=0; j<Ny; j++)
            {
                out_file << i <<", "<< j <<", "<<ux[j][i]<<", "<<uy[j][i]<<", "<<rho[j][i]<<", 1\n";
            }
        }
    out_file.close();
}