
/****** OE5450-Lid Driven Cavity Flow-Simple Algorithm Approach ******/
/****** OE21M025_Durga Rao ********/
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main()
{
    double u[55][55],u_st[55][55],d_e[55][55],v[55][55],v_st[55][55],d_n[55][55],p[55][55],p_st[55][55],p_c[55][55],
            b[55][55],u_n[55][55],v_n[55][55],p_n[55][55],v_f[55][55],u_f[55][55],p_f[55][55];
    double dx,dy,ae,aw,an,as,ap,ue,uw,vn,vs,error,alpha,dt,Re,aee,Aee,ann,Ann;

     int i,j,step;
     step=1;
     alpha=0.8;
     dx=0.02;
     dy=0.02;
     Re=100;
     error=1.0;


     //initializing variables
     //Final co-located grid
     for(i=1;i<52;i++)
        {
        for(j=1;j<52;j++){
            u_f[i][j]=0;
            v_f[i][j]=0;
            p_f[i][j]=1;

            u_f[1][j]=1;
        }
     }

     //staggered grid
     for(i=1;i<=52;i++)
        {
        for(j=1;j<52;j++){
                u[i][j]=0;
                u[1][j]=2;
                u_st[i][j]=0;
                d_e[i][j]=0;

                u_n[i][j]=0;
                u_n[1][j]=2;

        }
     }

     for(i=1;i<52;i++)
        {
        for(j=1;j<=52;j++){
                v[i][j]=0;
                 v_st[i][j]=0;
                 v_n[i][j]=0;
                 d_n[i][j]=0;
          }
        }

        for(i=1;i<=2;i++)
        {
             for(j=1;j<=52;j++){
                p[i][j] = 1;
                p_st[i][j] = 1;
                p_c[i][j]=0;
                p_n[i][j]=1;
                b[i][j]=0;

             }
        }

      ofstream eStream,sStream;
      eStream.open("output_err.txt");
      sStream.open("output_s.txt");

    while(error>1e-7)
    {
                    for(i=2;i<52;i++)
                        {
                        for(j=2;j<51;j++){
                            ue=0.5*(u[i][j]+u[i][j+1]);
                            uw=0.5*(u[i][j]+u[i][j-1]);
                            vn=0.5*(u[i-1][j]+v[i-1][j+1]);
                            vs=0.5*(u[i][j]+v[i][j+1]);
                         }
                        }

                   for(i=2;i<52;i++)
                        {
                           for(j=2;j<51;j++){


                            ae=-0.5*ue*dx+(1/Re);
                            aw=0.5*uw*dx+(1/Re);
                            an=-0.5*vn*dx+(1/Re);
                            as=0.5*vs*dx+(1/Re);
                           }
                        }
                    for(i=2;i<52;i++)
                        {
                           for(j=2;j<51;j++){
                            aee=0.5*dx*(ue-uw-vn-vs)+(4/Re);
                            Aee=-dy;
                            d_e[i][j]=(Aee/aee);

                            u_st[i][j]=(ae*u[i][j+1]+aw*u[i][j-1]+an*u[i-1][j]+as*u[i+1][j])/aee+d_e[i][j]*(p[i][j+1]-p[i][j]);


                    }
                }
                //x-momentum BC's
                     for(i=1;i<52;i++){
                       u_st[1][i]=2-u_st[2][i];
                       u_st[52][i]=-u_st[51][i];
                 }
                     for(i=2;i<52;i++){
                        u_st[i][1]=0;
                        u_st[i][51]=0;
                }

       //y-momentum equation interior
           for(i=2;i<51;i++)
            {
                for(j=2;j<52;j++){
                    ue=0.5*(u[i][j]+u[i+1][j]);
                    uw=0.5*(u[i][j-1]+u[i+1][j-1]);
                    vn=0.5*(u[i-1][j]+v[i][j]);
                    vs=0.5*(u[i][j]+v[i+1][j]);
                }
            }

            for(i=2;i<51;i++)
             {
                for(j=2;j<52;j++){

                    ae=-0.5*ue*dx+(1/Re);
                    aw=0.5*uw*dx+(1/Re);
                    an=-0.5*vn*dy+(1/Re);
                    as=0.5*vs*dy+(1/Re);
                }
             }

             for(i=2;i<51;i++)
               {
                for(j=2;j<52;j++){

                    ann=0.5*dx*(ue-uw-vn-vs)+(4/Re);

                    Ann=-dx;
                    d_n[i][j]=Ann/ann;

                    v_st[i][j]=(ae*v[i][j+1]+aw*v[i][j-1]+an*v[i-1][j]+as*v[i+1][j])/ann+d_n[i][j]*(p[i][j]-p[i+1][j]);

                }

           }

           //y-momentum BC'S
           for(i=1;i<52;i++){
            v_st[i][1]=-v_st[i][2];
            v_st[i][52]=-v_st[i][51];
           }
           for(j=2;j<52;j++){
            v_st[1][j]=0;
            v_st[51][j]=0;
           }

           //zeroing the corrections
           for(i=1;i<52;i++)
            {
                for(j=1;j<52;j++){
                        p_c[i][j]=0;
                }
           }

           //continuity equation aka pressure correction interior

           for(i=2;i<52;i++)
            {
                for(j=2;j<52;j++){
                        ae=-d_e[i][j]*dx;
                        aw=-d_e[i][j-1]*dx;
                        an=-d_n[i-1][j]*dy;
                        as=-d_n[i][j]*dy;

                        ap=ae+aw+an+as;

                        b[i][j]=-(u_st[i][j]-u[i][j-1])*dx+(v_st[i][j]-v_st[i-1][j])*dy;
                        p_c[i][j]=(ae*p_c[i][j+1] + aw*p_c[i][j-1] + an*p_c[i-1][j]+as*p_c[i+1][j]+b[i][j])/ap;
                }
           }

           //correction of the pressure field

               for(i=2;i<52;i++)
                {
                  for(j=2;j<52;j++){
                        p_n[i][j] = p[i][j] + alpha*p_c[i][j];
                  }
               }

          //continuity equation BC'S
              for(j=1;j<52;j++)
                {
                p_n[1][j] = p_n[2][j];
                p_n[52][j] = p_n[51][j];
              }
              for(i=1;i<52;i++){
                p_n[i][1] = p_n[i][2];
                p_n[i][52] = p_n[i][51];
              }

              //correcting velocities
              for(i=2;i<52;i++)
                {
                        for(j=2;j<51;j++){
                                u_n[i][j] = u_st[i][j] + alpha*d_e[i][j]*(p_c[i][j+1] - p_c[i][j]);
                        }
              }
               //x-momentum BC's
                     for(i=1;i<52;i++){
                       u_n[1][i]=2-u_n[2][i];
                       u_n[52][i]=-u_n[51][i];
                 }
                     for(i=2;i<52;i++){
                        u_n[i][1]=0;
                        u_n[i][51]=0;
                }

                for(i=2;i<51;i++){
                      for(j=2;j<52;j++){
                            v_n[i][j]=v_st[i][j]+alpha*d_n[i][j]*(p_c[i][j]-p_c[i+1][j]);
                      }
                }

          //y-momentum BC'S
           for(i=1;i<52;i++){
            v_n[i][1] = -v_n[i][2];
            v_n[i][52] = -v_n[i][51];
           }
           for(j=2;j<52;j++){
            v_n[1][j] = 0;
            v_n[51][j] = 0;
           }

           //continuity residual error measurement
           error = 0.0;
            for(i=2;i<51;i++)
                {
                      for(j=2;j<52;j++){
                            b[i][j] = -(u_st[i][j]-u[i][j-1])*dx+(v_st[i][j]-v_st[i-1][j])*dy;
                            error = error + fabs(b[i][j]);
                      }
            }

            if(step%1000==1)
                {
                cout<<"Error is "<<error<<" for the step "<<step<<endl;
                eStream<<error<<endl;
                sStream<<step<<endl;

               }
               eStream.close();
              sStream.close();

            //iterations

            for(i=2;i<52;i++)
                {
                        for(j=2;j<51;j++){
                                u[i][j]=u_n[i][j];
                        }
            }

            for(i=2;i<51;i++)
                {
                      for(j=2;j<52;j++){
                            v[i][j]=v_n[i][j];
                      }
            }

            for(i=2;i<52;i++)
                {
                  for(j=2;j<52;j++){
                        p[i][j]=p_n[i][j];
                  }
            }

            step = step + 1;

    }

      ofstream oStream;

    // After the converged solution, we map the staggered variables to collocated variables
    //u-velocity
     oStream.open("u.txt");
     for ( i = 2; i < 52; i++)
    {
        for ( j = 2;  j < 52; j++){
                u_f[i][j] = 0.5 * (u[i][j] + u[i+1][j]);
                oStream<<u_f[i][j]<<" ";
        }
        oStream<<endl;
    }
    oStream.close();
    cout<<"\nu_velocity output stored in path"<<endl;

    //v_velocity
    oStream.open("v.txt");
    for ( i = 1; i < 52; i++)
    {
        for ( j = 1;  j < 52; j++){
                v_f[i][j]=0.5 * (u[i][j]+u[i][j+1]);
                oStream<<v_f[i][j]<<" ";
        }
        oStream<<endl;
    }
    oStream.close();

    cout<<"\nv_velocity output stored in path"<<endl;

    //pressure-p
    oStream.open("p.txt");
    for ( i = 1; i < 52; i++)
    {
        for ( j = 1;  j < 52; j++){
            p_f[i][j] = 0.25 * (p[i][j] + p[i][j+1] + p[i+1][j] + p[i+1][j+1]);
           oStream<<p_f[i][j]<<" ";
        }
        oStream<<endl;
    }
    oStream.close();
    cout<<"\npressure output stored in path"<<endl;

    //centerline velocity
    oStream.open("center_u.txt");
    for(i=1; i<52; i++)
    {
        oStream<<u_f[i][26]<<endl;
    }
    oStream.close();

    cout<<"\ncenterline u_velocity output stored in path"<<endl;

    cout<<"\n-by OE21M025"<<endl;
  return 0;

}
