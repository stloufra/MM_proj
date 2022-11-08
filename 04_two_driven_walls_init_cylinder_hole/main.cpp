
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
const int Nx1 = Nx +1 ;
const int Ny1 = Ny +1 ;
const int L = Ny +1 ;
const int Nvel = 9 ; // Q number of possible descrete velocities
const double rho0 = 1.0;
const double ux0 = 0.0;
const double uy0 = 0.0;
const double uw = 0.1;
const double Re = 400.0;

//Model parameters

const int cx_pos[Nvel] = { 0, 1, 0, -1, 0, 1, -1, -1,  1};
const int cy_pos[Nvel] = { 0, 0, 1, 0, -1, 1,  1, -1, -1};
const double weight[Nvel] ={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};

//Cylinder
const double D_cylinder = 14 ;
const double x_cylinder = 100 ;
const double y_cylinder = 50 ;
bool objects [Ny1][Nx1];

//declaration of array of distribution function

double f[Ny1][Nx1][Nvel]; 
double f_post_col[Ny1][Nx1][Nvel]; 

//declaration of arrays of variables

double rho[Ny1][Nx1]; 
double ux[Ny1][Nx1]; 
double uy[Ny1][Nx1]; 

double ux0er[Ny1][Nx1];
double uy0er[Ny1][Nx1];

//help

double tau;
const int plot_every = 25;

//functions 

void initialization(void) ;
double f_equilibrium(double rho_eq, double ux_eq, double uy_eq, int k_eq);
void collision(void);
void streaming(void);
void postpro(void);
void bounce_back(void);
double Err(void);
void output(int s);

//help functions

double distance(int x1, int y1,int x2, int y2);
void boundary_objects(void);
void bo_bounce_back(void);
void bo_postpro(void);


// -----------------------------------------------------------------------------------------
// ---------------------------------------- MAIN -------------------------------------------
// -----------------------------------------------------------------------------------------

int main()
{
    int k,X2, Y2;
    double err;
    X2=Nx/2; Y2=Ny/2;
    k=0;
    err=1.0;
    tau = 3*L*uw/Re +0.5;

    initialization(); 

    while(k<5000)
    {
        k++;
        collision();
        streaming();
        bounce_back();
        postpro();

        if(k%1000==0)
        {
            err=Err();
            printf("err=%e ux_center=%e uy_center=%e rho_center=%e k=%d\n",err,ux[Y2][X2],uy[Y2][X2],rho[Y2][X2], k);
        }

        if(k%plot_every==0)
        {
            output(k);
        }

    }
    output(k);
    return 0;
}    

// -----------------------------------------------------------------------------------------
// -------------------------------- HELPING FUNCTIONS --------------------------------------
// -----------------------------------------------------------------------------------------

void boundary_objects()
{
    //structured bolean mesh values of existencion of objects
    int i, j, k;
    for (j=0; j<=Ny; j++)
    {
        for (i=0; i<=Nx; i++)
        {   
            if(distance(x_cylinder,y_cylinder,i,j) < D_cylinder)
            {
                objects[j][i]=true;
            }
            else
            {
                objects[j][i]=false;
            }
        }
     }
}

double distance(int x1, int y1,int x2, int y2)
{
    //simple distance calculation
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

void bo_bounce_back()
{
    int i, j;
 //boundary objects
    for (i=0; i<=Nx; i++)
    {
        for (j=0; j<=Ny; j++)
        { 
            if(objects[j][i])
            {
                f[j][i][3]=f_post_col[j][i][1];
                f[j][i][4]=f_post_col[j][i][2];
                f[j][i][1]=f_post_col[j][i][3];
                f[j][i][2]=f_post_col[j][i][4];
                f[j][i][7]=f_post_col[j][i][5];
                f[j][i][8]=f_post_col[j][i][6];
                f[j][i][5]=f_post_col[j][i][7];
                f[j][i][6]=f_post_col[j][i][8];

                
            }          
        } 
    }
}

void bo_postpro()
{
    int i, j;
for (j=0; j<=Ny; j++)
    {
        for (i=0; i<=Nx; i++)
        { 
            if(objects[j][i])
            {
                ux[j][i] = 0.0;
                uy[j][i] = 0.0;
                rho[j][i]= 0.0;
            }
        
        }

    }

}
// -----------------------------------------------------------------------------------------
// ----------------------------- EVOLUTION FUNCTIONS  --------------------------------------
// -----------------------------------------------------------------------------------------

void initialization() 
{
    // initialization with equilibrium method
    int j, i, k;

    boundary_objects();

    for (j=0; j<=Ny; j++)
    {
        for (i=0; i<=Nx; i++)
        {
            rho [j][i] = rho0;
            ux [j][i] = ux0;
            uy [j][i] = uy0;

            for(k=0;k<Nvel;k++)
            {
                f[j][i][k] = f_equilibrium(rho [j][i], ux[j][i], uy[j][i], k);
                
                //boundary of object
                if(objects[j][i])
                {
                    f[j][i][k] = 0;
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
    

    for (j=0; j<=Ny; j++)
    {
        for (i=0; i<=Nx; i++)
        {    
             for(k=0;k<Nvel;k++)
            {
                feq = f_equilibrium(rho [j][i], ux[j][i], uy[j][i], k);
                f_post_col[j][i][k]= f[j][i][k] - (f[j][i][k]-feq)/tau;                             //collision
            }
        }
    }
}

void streaming()
{
    int j, i, jd, id, k;
    for (j=0; j<=Ny; j++)
    {
        for (i=0; i<=Nx; i++)
        {    
             for(k=0;k<Nvel;k++)
            {
                jd = j - cy_pos[k];
                id = i - cx_pos[k];
                
                if(jd>=0 && jd<=Ny && id>=0 && id<=Nx)
                {
                    f[j][i][k]=f_post_col[jd][id][k];                                               //streaming
                }

            }
        }
    }
}

void bounce_back()
{
    int i,j;
    
    //upper wall
    for(i=0;i<=Nx; i++)
    {
        f[Ny][i][4]=f_post_col[Ny][i][2];
        f[Ny][i][7]=f_post_col[Ny][i][5]+6*rho[Ny][i]*weight[7]*cx_pos[7]*uw;
        f[Ny][i][8]=f_post_col[Ny][i][6]+6*rho[Ny][i]*weight[8]*cx_pos[8]*uw;
     }
    
    //lower wall
    for(i=0;i<=Nx;i++)
    {
        f[0][i][2]=f_post_col[0][i][4];
        f[0][i][5]=f_post_col[0][i][7]+6*rho[0][i]*weight[5]*cx_pos[5]*uw;
        f[0][i][6]=f_post_col[0][i][8]+6*rho[0][i]*weight[6]*cx_pos[6]*uw;
    }
    //left wall
    for(j=0;j<=Ny;j++)
    {
        f[j][0][1]=f_post_col[j][0][3];
        f[j][0][5]=f_post_col[j][0][7];
        f[j][0][8]=f_post_col[j][0][6];
    }

    //right wall
    // i=Nx: right wall
    for(j=0;j<=Ny;j++)
    {
        f[j][Nx][3]=f_post_col[j][Nx][1];
        f[j][Nx][7]=f_post_col[j][Nx][5];
        f[j][Nx][6]=f_post_col[j][Nx][8];
    }
    
    //bo_bounce_back();

}
void postpro()
{
    int i, j;
    for (j=0; j<=Ny; j++)
    {
        for (i=0; i<=Nx; i++)
        { 
            rho[j][i]=f[j][i][0]+f[j][i][1]+f[j][i][2]+f[j][i][3]+f[j][i][4]+f[j][i][5]+f[j][i][6]+f[j][i][7]+f[j][i][8];
            ux[j][i]=(f[j][i][1]+f[j][i][5]+f[j][i][8]-f[j][i][3]-f[j][i][6]-f[j][i][7])/rho[j][i];
            uy[j][i]=(f[j][i][5]+f[j][i][6]+f[j][i][2]-f[j][i][7]-f[j][i][8]-f[j][i][4])/rho[j][i];      
        }

    }
    //bo_postpro();
}

double Err() 
{
    int j, i;
    double er1, er2;
    er1=0.0; er2=0.0;

    for (j=1; j<Ny; j++)
    {
        for (i=0; i<Nx; i++)
        {
            er1+=sqrt((ux[j][i]-ux0er[j][i])*(ux[j][i]-ux0er[j][i])+(uy[j][i]-uy0er[j][i])*(uy[j][i]-uy0er[j][i]));
            er2+=sqrt(ux[j][i]*ux[j][i]+uy[j][i]*uy[j][i]);
            ux0er[j][i]=ux[j][i];
            uy0er[j][i]=uy[j][i];

           
        }
    }
    return er1/er2;
}

void output(int s)
{
    int i, j;

    double u[Ny1][Nx1];
    double curl[Ny1][Nx1];

    for (i=0; i<=Nx; i++)
    {
        for (j=0; j<=Ny; j++)
        {
            u[j][i]=sqrt(ux[j][i]*ux[j][i]+uy[j][i]*uy[j][i]);
        }
    }

    
    std::string step =std::to_string(s/plot_every);

    std::ofstream out_file("LBE.csv."+step);

    out_file << "x coord, y coord, ux, uy, u, rho, obj\n" ;
    
    for (i=0; i<Nx; i++)
        {
            for (j=0; j<Ny; j++)
            {
                out_file << i <<", "<< j <<", "<<ux[j][i]<<", "<<uy[j][i]<<", "<<u[j][i]<<", "<<rho[j][i]<<", "<<objects[j][i]<<", 1\n";
            }
        }

    out_file.close();


}

