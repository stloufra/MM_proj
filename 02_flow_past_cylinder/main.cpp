#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

// -----------------------------------------------------------------------------------------
// ------------------------ CONSTANTS AND DECLARATION --------------------------------------
// -----------------------------------------------------------------------------------------

const int Nx = 400 ;
const int Ny = 100 ;
const int Nx1 = Nx +1 ;
const int Ny1 = Ny +1 ;

const int Nt = 300 ;
const int plot_every = 1;
double tau;

// cylinder
const double D_cylinder = 14 ;
const double x_cylinder = 100 ;
const double y_cylinder = 50 ;
bool objects [Ny1][Nx1];


//using D2Q9 model

const int Nvel = 9 ; // number of possible descrete velocities

const int cx_pos[Nvel] = { 0, 1, 0, -1, 0, 1, -1, -1,  1};
const int cy_pos[Nvel] = { 0, 0, 1, 0, -1, 1,  1, -1, -1};
const double weight[Nvel] ={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};

//initial conditions

const double rho0 = 1.0;
const double ux0 = 1.0;
const double uy0 = 0.0;

double const ny = 1.035e-6; //water

//declaration of arrays of variables

double rho[Ny1][Nx1]; 
double ux[Ny1][Nx1]; 
double uy[Ny1][Nx1]; 

//declaration of array of distribution function

double f[Ny1][Nx1][Nvel]; 
double f_post_col[Ny1][Nx1][Nvel]; 

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
                f[j][i][1] = 3.4;

                //boundary of object
                if(objects[j][i])
                {
                    f[j][i][k] = 0;
                }
            }

        }
    }    

return 0;
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

int bounce_back()
{
    int i,j;

   //lower wall
    for(i=0;i<=Nx;i++)
    {
        f[0][i][2]=f_post_col[0][i][4];
        f[0][i][5]=f_post_col[0][i][7];
        f[0][i][6]=f_post_col[0][i][8];
    }
    //left wall
    for(j=0;j<=Ny;j++)
    {
        f[j][0][1]=f_post_col[j][0][3];
        f[j][0][5]=f_post_col[j][0][7];
        f[j][0][8]=f_post_col[j][0][6];
    }

    //boundary objects
    for (i=0; i<Nx; i++)
    {
        for (j=0; j<Ny; j++)
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

                //inside cylinder 
            }
        }
}

void test(int t)
{
    std::cout << "\n f > "<< f[20][200][1] << "\n";
    std::cout << "\n rho > "<< rho[20][200] << "\n";
    std::cout << "\n ux > "<< ux[20][200] << "\n";

    std::cout << "Time > "<< t << "\n";
}

void output(int s)
{
    int i, j;

    std::string step =std::to_string(s);

    std::ofstream out_file("LBE.csv."+step);

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

// -----------------------------------------------------------------------------------------
// ---------------------------------------- MAIN -------------------------------------------
// -----------------------------------------------------------------------------------------

int main()
{
    int t,s=0;
    tau = tau_f();
    initialization_by_u();

    for (t=0;t<Nt;t++)
    {
        collision();
        streaming();
        bounce_back();
        postpro();
        test(t);
        if(t%plot_every==0)
        {
            output(s);
            s++;
        }


    }
    return 0;
}
