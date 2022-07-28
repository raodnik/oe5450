/***** Lagrange Interpolation *****/

#include<bits/stdc++.h>

using namespace std;

/**** Function declaration  *****/

double L_interpolation_V(double x[], double y[], double m);
double L_interpolation_A(double x[], double y[], double m);

  /**** Function call  *****/

int main()
{
   double x[21] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};      //data set
   double y[21] = {0,2,7,12,22,35,48,53,67,90,120,87,62,35,23,18,15,9,3,1,0};
   double c, a=1, b=3;
   c=a/b;
   while(c<=20)
   {
       cout<<"V("<<c<<"):"<<L_interpolation_V(x, y, c)<<"   "<<"A("<<c<<"):"<<L_interpolation_A(x, y, c)<<endl;
       c++;
   }
   cout<<"\n-by OE21M025"<<endl;

}

 /***** Function definition for Velocity *****/

double L_interpolation_V(double x[], double y[], double m)
{
    double s=0;
    for(int i=0; i<=20; i++)
    {
        double p=1;
        for(int j=0; j<=20; j++)
        {
            if(i!=j)
            {
            p *= (m-x[j])/(x[i]-x[j]);
            }
        }
        p *= y[i];
        s += p;

    }
    return s;
}

   /***** Function definition for Acceleration *****/

double L_interpolation_A(double x[], double y[], double m)
{
    double q=0,s=0;
    for(int i=0; i<=20; i++)
    {
        for(int j=0; j<=20; j++)
        {
            double k=1;
                for(int g=0; g<=20; g++)
                {
                    if(g!=i && g!=j){
                        k*=(m-x[g])/(x[i]-x[g]);
                    }
                }
                if(i!=j){
                q+=k/(x[j]-x[i]);
            }

        }
         s += y[i]*q;

    }
    return s;
}

