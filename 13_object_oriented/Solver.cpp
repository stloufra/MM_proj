#include "Solver.h"

Solver::Solver(int Ny_in, int Nx_in, Mesher mesh_in):
rho(Ny_in,Nx_in),
ux(Ny_in,Nx_in),
uy(Ny_in,Nx_in),
df(Ny_in,Nx_in),
df_post(Ny_in,Nx_in),
mesh(Ny_in,Nx_in),
ux0er(Ny_in,Nx_in),
uy0er(Ny_in,Nx_in)
{
    mesh = mesh_in;
    Nx = Nx_in;
    Ny = Ny_in;
}

//--------------------------INITIALIYATION---------------------------------------


void Solver::initialization_eq(double rho0, double ux0, double uy0,double Re)
{
    for(int j=0; j < Ny; j++)
  {
    for(int i=0; i < Nx; i++)
    {
        if(mesh.objects(j,i)==1)
        {
        rho(j,i)=rho0;
        ux(j,i) = ux0;
        uy(j,i) = uy0;
        for(int k=0;k<Nvel;k++)
        {
            df(j,i,k)=f_equilibrium(rho(j,i), ux(j,i), uy(j,i), k);
        }
        }

        if(mesh.objects(j,i)==0)
        {
        rho(j,i)=0;
        ux(j,i)=0;
        uy(j,i)=0;
        for(int k=0;k<Nvel;k++)
        {
            df(j,i,k)=0;
        }
        }
    }
  }
    tau = 3*Nx*std::max(sqrt(ux0*ux0+uy0*uy0),0.1)/Re +0.5;

    assert(tau>0.51 && tau<1);

    std::cout<<"Initialization undergone successfully."<<"\n";

}

//--------------------------EVOLUTION FUNCTIONS---------------------------------------

void Solver::collision() 
{
    
    double feq;
    

    for (int j=0; j<Ny; j++)
    {
        for (int i=0; i<Nx; i++)
        {   
            if(mesh.objects(j,i)==1)
            {
                for(int k=0;k<Nvel;k++)
                {
                    feq = f_equilibrium(rho(j,i), ux(j,i), uy(j,i), k);
                    df_post(j,i,k)= df(j,i,k) - (df(j,i,k)-feq)/tau;                             //collision
                }
            }
        }
    }
}

void Solver::streaming()
{
    int j, i, jd, id, k;
    for (j=0; j<Ny ; j++)
    {
        for (i=0; i<Nx; i++)
        {    
            for(k=0;k<Nvel;k++)
            {
                if(mesh.objects(j,i)==1)
                {
                    jd = j - cy_pos[k];
                    id = i - cx_pos[k];
                    if(jd>=0 && jd <= Ny && id>=0 && id <= Nx)
                    df(j,i,k)=df_post(jd,id,k);                                               //streaming
                    
                }

                if(mesh.objects(j,i)==0)
                {
                   df(j,i,k)=0;
                } 
            }
            
        }
    }
}

void Solver::bounce_back()
{
    
    for (int j=0; j<Ny; j++)
    {
        for (int i=0; i<Nx; i++)
        {
            //rovne stacionarni steny
            if(mesh.mesh(j,i)==2)
            {
                df(j,i,2)=df_post(j,i,4);
                df(j,i,5)=df_post(j,i,7);
                df(j,i,6)=df_post(j,i,8);
            }

            if(mesh.mesh(j,i)==3)
            {
                df(j,i,3)=df_post(j,i,1);
                df(j,i,6)=df_post(j,i,8);
                df(j,i,7)=df_post(j,i,5);
            }

            if(mesh.mesh(j,i)==4)
            {
                df(j,i,4)=df_post(j,i,2);
                df(j,i,7)=df_post(j,i,5);
                df(j,i,8)=df_post(j,i,6);
            }

            if(mesh.mesh(j,i)==5)
            {
                df(j,i,1)=df_post(j,i,3);
                df(j,i,5)=df_post(j,i,7);
                df(j,i,8)=df_post(j,i,6);
            }

            //rohy vnitrni
            if(mesh.mesh(j,i)==6)
            {
                df(j,i,2)=df_post(j,i,4);
                df(j,i,1)=df_post(j,i,3);
                df(j,i,5)=df_post(j,i,7);
            }

            if(mesh.mesh(j,i)==7)
            {
                df(j,i,2)=df_post(j,i,4);
                df(j,i,3)=df_post(j,i,1);
                df(j,i,6)=df_post(j,i,8);
            }

            if(mesh.mesh(j,i)==8)
            {
                df(j,i,4)=df_post(j,i,2);
                df(j,i,3)=df_post(j,i,1);
                df(j,i,7)=df_post(j,i,5);
            }

            if(mesh.mesh(j,i)==9)
            {
                df(j,i,4)=df_post(j,i,2);
                df(j,i,1)=df_post(j,i,3);
                df(j,i,8)=df_post(j,i,6);
            }

            //rohy vnejsi
            if(mesh.mesh(j,i)==10)
            {
                df(j,i,7)=df_post(j,i,5);
            }

            if(mesh.mesh(j,i)==11)
            {
                df(j,i,8)=df_post(j,i,6);
            }

            if(mesh.mesh(j,i)==12)
            {
                df(j,i,5)=df_post(j,i,7);
                
            }

            if(mesh.mesh(j,i)==13)
            {
                df(j,i,6)=df_post(j,i,8);
            }

            //rovne pohyblive steny
            if(mesh.mesh(j,i)==22)
            {
                df(j,i,2)=df_post(j,i,4);
                df(j,i,5)=df_post(j,i,7)+6*rho(j,i)*weight[5]*cx_pos[5]*mesh.uw;
                df(j,i,6)=df_post(j,i,8)+6*rho(j,i)*weight[6]*cx_pos[6]*uw22; 
            }

            if(mesh.mesh(j,i)==24)
            {
                df(j,i,4)=df_post(j,i,2);
                df(j,i,7)=df_post(j,i,5)+6*rho(j,i)*weight[7]*cx_pos[7]*uw24;
                df(j,i,8)=df_post(j,i,6)+6*rho(j,i)*weight[8]*cx_pos[8]*uw24;
            }
        }
    }

   
    
}

void Solver::postpro()
{
    int i, j;
    for (j=0; j<Ny; j++)
    {
        for (i=0; i<Nx; i++)
        { 
            if(mesh.objects(j,i)==1)
            {
                rho(j,i)=df(j,i,0)+df(j,i,1)+df(j,i,2)+df(j,i,3)+df(j,i,4)+df(j,i,5)+df(j,i,6)+df(j,i,7)+df(j,i,8);
                ux(j,i)=(df(j,i,1)+df(j,i,5)+df(j,i,8)-df(j,i,3)-df(j,i,6)-df(j,i,7))/rho(j,i);
                uy(j,i)=(df(j,i,5)+df(j,i,6)+df(j,i,2)-df(j,i,7)-df(j,i,8)-df(j,i,4))/rho(j,i);      
            }
        }

    }
}

void Solver::Err() 
{
    int j, i;
    double er1, er2;
    er1=0.0; er2=0.0;

    for (j=1; j<Ny-1; j++)
    {
        for (i=0; i<Nx-1; i++)
        {
            er1+=sqrt((ux(j,i)-ux0er(j,i))*(ux(j,i)-ux0er(j,i))+(uy(j,i)-uy0er(j,i))*(uy(j,i)-uy0er(j,i)));
            er2+=sqrt(ux(j,i)*ux(j,i)+uy(j,i)*uy(j,i));
            ux0er(j,i)=ux(j,i);
            uy0er(j,i)=uy(j,i);

           
        }
    }
    err = er1/er2;
}
//--------------------------HELP FUNCTIONS---------------------------------------

double Solver::f_equilibrium(double rho_eq, double ux_eq, double uy_eq, int k_eq)
{
    double cu, u2;

    cu = cx_pos[k_eq]*ux_eq + cy_pos[k_eq]*uy_eq;
    u2 = ux_eq*ux_eq + uy_eq*uy_eq;

    return weight[k_eq]*rho_eq*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u2);                                    //equlibrium
}

void Solver::output_VTK(int s,int plot_every)
{
       
    std::string step =std::to_string(s/plot_every);

    std::ofstream out_file("results/LBE."+step+".vtk");

    out_file << "# vtk DataFile Version 2.0\n" ;
    out_file << "LBE mesh\n" ;
    out_file << "ASCII\n";
    out_file << "DATASET STRUCTURED_POINTS\n" ;
    out_file << "DIMENSIONS "<<Nx<<" "<<Ny<<" 1\n" ;
    out_file << "ASPECT_RATIO 1 1 1\n";
    out_file << "ORIGIN 0 0 0\n";
    out_file << "POINT_DATA "<<Nx*Ny<<"\n";

    out_file << "SCALARS "<<"rho "<<"double 1\n";
    out_file << "LOOKUP_TABLE default\n";
    for(int j=0; j < Ny; j++)
    {
        for(int i=0; i < Nx; i++)
        {
            out_file <<rho(j,i)<<"\n";
        }
    }
       
    out_file << "VECTORS "<<"U "<<"double\n";
    for(int j=0; j < Ny; j++)
    {
        for(int i=0; i <Nx; i++)
        {
            out_file <<ux(j,i)<<" "<<uy(j,i)<<" 0\n";
        }
    }


    out_file.close();


}

Solver::~Solver()
{

}