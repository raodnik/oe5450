/* C++ CODE for Computing the Area of triangle using the approximation
of triangle into squares and rectangles and find the
accuracy */


#include <bits/stdc++.h>
#include<cmath>

using namespace std;

int main()
{
   float b,h,x;                          // Declaring input Variables
   double m,n,s,t,A1,A2,A3,A4,A5,K,L;   //Declaring the variables to be estimated
   cout<<"Enter the base width: "<<endl;
   cin>>b;
   cout<<"Enter the height: "<<endl;
   cin>>h;

   n=b/1.5;                 // base width for first triangle
   m=0.5*n;                 //base width of second triangle

   s=sqrt(pow(m,2)+pow(h,2));
   t=sqrt(pow(n,2)+pow(h,2));

   A1=0.5*(m*h+n*h);
   cout<<"Area of triangle by analytical method: "<<A1<<endl;

   float i;
   x=1;                       //initializing the grid number
   cout<<"Enter number of grids: "<<endl;
   cin>>x;

   for(i=1; i<x; i++)
   {
       A2 += (h*(1-i/x)*(n/x));
       A4 += (h*(1-i/x)*(m/x));
   }

   A3 = (n/(2*x));
   A5 = (m/(2*x));
   cout<<"Area of type-1 rectangle in first triangle: "<<A2<<endl;
   cout<<"Area of type-1 rectangle in second triangle: "<<A4<<endl;
   cout<<"Area of type-2 rectangle in first triangle: "<<A3<<endl;
   cout<<"Area of type-2 rectangle in second triangle: "<<A5<<"\n"<<endl;

   K=A2+A3+A4+A5;                // Area of triangle using numerical method
   cout<<"Area of triangle using numerical method: "<<K<<"\n"<<endl;

   L=((A1-K)/A1)*100;           //computing Accuracy
   cout<<"Error % is: "<<L<<endl;

}


/* END OF THE PROGRAM */
