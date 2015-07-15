//
//  sim8.c
//  
//
//  Created by Arthur Romero on 7/10/15.
//
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <complex>
#include <ccomplex>
#include "complex.h"


#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TROOT.h"


#include "SimplePlot.hh"


//global variables:


double e_charge = 0.00000000000000000016021765;
double e_mass = 0.00000000000000000000000000000091093819;
double Plancks = 0.0000000000000000000000000000000006626068;
double Coulomb = 0.00000000898755179;
double pi = 3.14159265358979;
double const_e = 2.71828182845905;
double chargeratio =0;//strength of image charge(units of e (?))

double lambda = 0.00000000001;

double eta1 = .4;//G1 open fraction
double eta2 = .4;//G2 open fraction
double d = 0.0000001;// period of grating

double r0 = -4.04;//initial radius of wavefront curvature
double el0 = 0.000001;// initial coherence width
double w0 = 0.00003;// initial beam width

double G1_z = 0.000001;
double G2_z = 1;
double G2_x = 0.00000005; //d/2;

double theta = 0;

int r =41;// length of the array.

double energy = 15000;
double width = 0.00000004;
double thick = 0.000000014;
double wedgeangle = 0;
double tilt =0;

double res = 1000;

double zstart = -0.1;
double zend = 2.1;
double xstart = -0.00020;
double xend = 0.00020;
double ystart = -0.00011;
double yend = 0.00011;

double xpnts = 300;
double ypnts = 300;
double zpnts = 300;
//prototype functions:




double zp(double z, double v);// prototype

double w(double z,double r0, double el0, double w0);// prototype

double el(double z, double r0, double el0, double w0);// prototype

double v(double z,double r0, double el0, double w0);// prototype

double sinc(double x);//prototype


void gp1(double z,double r0,double el0, double w0);//prototype

void gp2(double z12,double z23, double mytheta, double el1x, double w1x, double r1x, double el1y, double w1y, double r1y, double G2_x); // prototype

int x2pnts(float *arr, int r, int value);//prototype













int main( )
{
    double lam = pow((150)/(4000),1/2)*pow(10,-10);// double lam = 1.936E-10//
    
    ///////////Fourier Components Variables//////////
    
    double res = 1000;// double res = 1000;
    double eta = width/d;
    double vel = pow(2*energy*e_charge/e_mass,1/2);
    double alpha = wedgeangle*pi/180;
    double beta = tilt*pi;// since tilt=0, beta=0.
    
    
    
    
    
    long double xmin;
    long double xmax;
    
    
    
    double zres = (zend-zstart)/zpnts;// matches with Mathematica code value
    //printf("the value of zres is: %f\n",zres);
    
    
    
    double w1;
    double r1;
    double el1;
    double el2;
    double *x1;
    

    
    
    w1=w(G1_z,r0,el0,w0);//gp1
    //printf("the value of w1 is:  %.12f\n",w1);// matches with Mathematica value of w1
    
    r1 = v(G1_z,r0,el0,w0);//gp1
    //printf("the value of r1 is:  %.12f\n",r1);//matches with Mathematica value of r1
    
    
    el1=el(G1_z,r0,el0,w0);//gp1
    //printf("the value of el1 is:  %.12f\n",el1); //matches with Mathematica value of el1
    
    el2 = el(zstart-G1_z+zres*100, r1, el1, w1);//gp2
     //printf("the value of el2 is:  %.12f\n",el2);// matches with mathematica value
    
    
    
     gp2(G2_z-G1_z,zstart+0*zres,theta,el1,w1,r1,el1,w1,r1,G2_x);
    
     gp1(zstart - G1_z + zres*100,r1,el1,w1);
    // printf("the value of el2 is:  %.12f\n",el2);
    
    
    
    

    
    
}






double zp(double z, double v)
{
    
    double zp;
    zp = (v*z)/(z+v);
    return(zp);
    
}



double w(double z,double r0, double el0, double w0)
{
    
    double w;
    double test4;
    double test5;
    double test6;
    test4 = (z)/(zp(z,r0));
    test5 = pow(((lambda*zp(z,r0))/(el0*w0)),2);
    test6 =(sqrt((1+test5)));
    
    w = (el0)*(fabs((z)/(zp(z,r0)))*((sqrt((1+test5)))));
    
    return(w);
    
}

double el(double z,double r0, double el0, double w0)
{
    
    double w;
    
    double test4;
    double test5;
    double test6;
    test4 = (z)/(zp(z,r0));
    test5 = pow(((lambda*zp(z,r0))/(el0*w0)),2);
    test6 =(sqrt((1+test5)));
    
    
    //w = el0*fabs((z)/(((zp(z,r0)))))*pow(1+pow(lambda,2)*pow(zp(z,r0),2)/((pow((el0*w0),2))),1/2);
    // w = el0*fabs((z)/(((zp(z,r0)))))*pow((1+pow((lambda*zp(z,r0)/(el0*w0)),2)),1/2);
    
    
    w = (el0)*(fabs(test4))*(test6);
    
    
    //printf("the value of test4 is : %0.12f\n",test4);
    //printf("the value of test5 is : %0.12f\n",test5);
    //printf("the value of test6 is : %0.12f\n",test6);
    
    
    
    
    
    
    
    return(w);
    
}

double v(double z,double r0, double el0, double w0)
{
    
    double v;
    
    
    v=(z)/(1-zp(z,r0)/(z*(1+pow(((lambda*zp(z,r0)/(el0*w0))),2))));
    
    return(v);
}



double sinc(double x)
{
    
    double sinc;
    
    if (x==0)
    {
        sinc = 1;
    }
    else
    {
        sinc = (sin(x))/x;
    }
    
    return(sinc);
}


int x2pnts(float *arr, int r, int value)// put number of rows and columns as arguments?
{
    
    int c = 2;
    int arr2 = {0};
    
    int ans;
    for (int i = 0; i < r; i++)
    {
        //printf("%f \n", *((arr+i*n) + j));
        arr2 = *(arr+i*c);
        //printf("the values of arr2 are: %d \t %d \n", arr2,i);
        
        if (arr2 == value) {
            //printf("the value of i is: %d\n",i+1);
            //printf("the equivalent for %d in ReT is the %d th element: %f\n",value,i+1, *((arr+(i)*2) + 1));
            ans = (*((arr+(i)*2) + 1));
            //printf("the value of ans is: %f\n", ans);
            return(i);
        }
    }
    
    
    
}

void gp1(double z,double r0,double el0, double w0)
{
    double coef;
    double cutoff=pow(10,-3);
    double lim=5;
    
    
    
    
    double eta = width/d;
    double vel = pow(2*energy*e_charge/e_mass,1/2);
    double alpha = wedgeangle*pi/180;
    double beta = tilt*pi;
    double w1;
    double r1;
    double el1;
    double el2;
    
    
    double zres = (zend-zstart)/zpnts;
    
    
    
    double w2;
    double r2;
    
    long double xmin;
    long double xmax;
    
    float fc;
    float ph;
    float ex;
    int j;
    
    
    
    w1=w(G1_z,r0,el0,w0);
    w2=w(z,r0,el0,w0);
    
    
    r1 = v(G1_z,r0,el0,w0);
    r2 = v(z,r0,el0,w0);
    //printf("the values of r2 and w2 are:  %f \t %f \n",r2,w2);
    
    
    el1=el(G1_z,r0,el0,w0);
    //printf(" the values of zstart,G1_z,zres, r1,el1 and w1 are: %f \t %f \t %f \t %f \t %f \t %f\n",zstart,G1_z,zres,r1,el1,w1);
    
    el2 = el(z, r1, el1, w1);
    //printf("the value of el2 is : %0.12f\n",el2);
    
    
    float ReT[41][2]={{0}};
    for (int i=0; i<=41; i++) {
        ReT[i][0]=i-20;
        
    }
    
    
    
    float ImT[41][2]={{0}};
    for (int i=0; i<=41; i++) {
        ImT[i][0]=i-20;
        
    }
    
    
    
    
    
    
    if (beta>=0){
        
        
        
        xmin= width*(1/res - cos(beta)/2);
        if (beta<=alpha) {
            xmax=(width*cos(beta))/2-width/res;
        }
        else
        {
            xmax=width*cos(beta)/2-width/res+thick*(tan(alpha)-tan(beta));
        }
    }
    else
    {
        xmax = (width*cos(beta)/2)-width/res;
        if (fabsl(beta)<=alpha) {
            xmin = -((width*cos(beta))/2)+width/res;
        }
        else
        {
            xmin = -((width*cos(beta))/2)+width/res - thick*(tan(alpha)-tan(beta));
        }
        
    }
    
    
    
    
    
    for(int n=-20;n<=20;n++)
    {
        for(ex=xmin; ex<xmax; ex+=width/res)//copied from above
        {
            fc = 2*pi*n*ex/d;
            
            ph = -width*thick*chargeratio*pow(e_charge,2)*(2*pi*Coulomb/Plancks)/(vel*(.25*pow(width,2)-pow(ex,2)));
            
            j=n+20;
            
            ReT[j][1] += cos(ph+fc);
            ImT[j][1] += sin(ph+fc);
            
        }
        
        
        
    }
    
    for (int i=0; i<=40; i++) {
        
        ReT[i][1] = ReT[i][1]/res;
        ImT[i][1] = ImT[i][1]/res;
        
    }
    
    
    
    
    double ix[300][2]={0};
    
    
    for (int i=0; i<300; i++) {
        ix[i][0]= xstart+(i-1)*((xend-xstart)/(xpnts-1));
        //printf("the values of ix[i] ddand i are: %f\t and %d\n",ix[i][1],i);//different format, but same values as in the Mathematica code
    }
    
    
    
    
    
    
    
    
    
    
    
    for (int i=0; i<300; i++) {
        for (int n=-lim; n<=lim; n++) {
            
            for (int m=-lim; m<=lim; m++) {
                double dn =n-m;
                double dm = (m+n)/2;
                double test1=0;
                double coef2=0;
                
                
                
                
                
                
                
                
                if (test1==1)
                {
                    coef = sinc(eta1*pi*n)*(sinc(eta1*pi*m)*pow((eta1), 2));
                    
                    //printf("the value of eta1 is: %f\n",eta1);
                }
                else
                {
                    coef =  ReT[x2pnts((float *)ReT, r, n)][1]*ReT[x2pnts((float *)ReT, r, m)][1]+ImT[x2pnts((float *)ImT, r, n)][1]*ImT[x2pnts((float *)ImT, r, m)][1];//
                    //printf("the value of eta1 is: %f\n",eta1);
                    //printf("the value %d \t and %d \t of ReT and ImT are: %0.9f \t %0.9f \n",(n+20),(m+20),ReT[(n+20)][1],ImT[(m+20)][1]);
                    
                    //printf("the value of coef, m and n are: %0.10f \t %d \t %d \n",coef,m,n);//!!!! results match with mathematica so far.
                    
                    
                }
                
                
                
                //printf("the value of el2 is:    %0.9f\n",el2);//ok, it matches.
                
                coef = coef*exp(-pi*(dn*(lambda*z/(pow(d*el2,2)))));// added isfinite macro in order to avoid inf values
                if (isfinite(coef)==0){
                    coef=0;
                }
                
                
                //printf("the value of coef and cutoff are: %0.10f \t %f \t %d \t %d \n",coef, cutoff,n,m,dn,dm);// coef values match
                
                
                
                if (coef>=cutoff) {
                    
                    ix[i][1] = ix[i][1] +  (coef*exp(-pi*pow(((ix[i][0]-dm*lambda*z/d)/w2),2)*cos(2*pi*(dn/d)*(ix[i][0]-dm*lambda*z/d)*(1-z/r2))));
                    
                    //printf("the value of coef2 is: %0.10f \t %d \t %d \t %f \t %f \n",coef2,n,m,dn,dm);//coef2 values are too low
                    
                    
                    //ix[i][1] = ix[i][1] + coef2;
                    
                    
                    
                    
                    continue;
                    
                }
            }
            
            
            
            
            
            
            
            
            
        }
        //printf("the values of ix[i][1] are:  %0.19f \t %d \n",ix[i][1],i);//
        
    }
    
    
    double ix1[300]={0};
    double ix2[300]={0};
    
    
    
    for (int i=0; i<300; i++) {
        for (int j=0; j<2; j++) {
            if (j==0) {
                ix1[i]=ix[i][j];
            }
            if (j==1) {
                ix2[i] = ix[i][j];
            }
        }
    }
    
    
    SimplePlot::graph("gp1 graph",ix1,ix2,300);
    
    
    
    
}











void gp2(double z12,double z23, double mytheta, double el1x, double w1x, double r1x, double el1y, double w1y, double r1y, double G2_x)
{
    
    
    
    
    double res = 1000;
    double eta = width/d;
    double vel = pow(2*energy*e_charge/e_mass,1/2);
    double alpha = wedgeangle*pi/180;
    double beta = tilt*pi;
    double w1;
    double r1;
    double el1;
    double el2;
    
    double theta = pi*mytheta/180;
    double d1=d;
    double d2=d;
    double z13 = z12+z23;
    
    long double xmin;
    long double xmax;
    
    float fc;
    float ph;
    float ex;
    
    
    double phi =0;
    
    double cutoff = 0.001;
    
    double lim =5;
    
    double _Complex coef;
    //double complex coef;
    double dn = 0;
    double dm =0;
    double m=0;
    double n=0;
    int a =0;
    int b =0;
    int c=0;
    int d5=0;
    double test1=0;
    

   
    
    
    
    // printf("the values of theta, d1,d2 and z13 are:   %0.15f \t %0.15f \t %0.15f \t %0.15f \n",theta,d1,d2,z13);// the values match with the mathematica values.
    
    
    
    double el3x = el(z13, r1x, el1x, w1x);//G2z - G1z + zstart + 0*zres, r1, el1, w1
    // printf("the value of el3x is: %0.15f\n",el3x);//value matches with mathematica
    
    
    double w3x = w(z13,r1x,el1x,w1x);
    //printf("the value of w3x is: %0.15f\n",w3x);//value matches with mathematica
    
    
    double v3x = v(z13,r1x,el1x,w1x);
    // printf("the value of v3x is: %0.15f\n",v3x);//value matches with mathematica
    
    
    double el3y = el(z13,r1y,el1y, w1y);
    // printf("the value of el3y is: %0.15f\n",el3y);//value matches with mathematica
    
    
    double w3y = w(z13,r1y,el1y,w1y);
    // printf("the value of w3y is: %0.15f\n",w3y);//value matches with mathematica
    
    
    double v3y = v(z13,r1y,el1y,w1y);
    //printf("the value of v3x is: %f\n",v3x);//value matches with mathematica
    
    
    
    
    
    
    
    float ReT[41][2]={{0}};
    for (int i=0; i<41; i++) {
        ReT[i][0]=i-20;
        //printf("value of ReT is: %f and j is: %d \n",ReT[i][1],(i-20)); ok
        
    }
    
    
    
    float ImT[41][2]={{0}};
    for (int i=0; i<41; i++) {
        ImT[i][0]=i-20;
        //printf("value of ImT is: %f and j is: %d \n",ImT[i][1],(i-20)); ok
        
    }
    
    
    
    if (beta>=0){
        
        
        
        xmin= width*(1/res - cos(beta)/2);
        if (beta<=alpha) {
            xmax=(width*cos(beta))/2-width/res;
        }
        else
        {
            xmax=width*cos(beta)/2-width/res+thick*(tan(alpha)-tan(beta));
        }
    }
    else
    {
        xmax = (width*cos(beta)/2)-width/res;
        if (fabsl(beta)<=alpha) {
            xmin = -((width*cos(beta))/2)+width/res;
        }
        else
        {
            xmin = -((width*cos(beta))/2)+width/res - thick*(tan(alpha)-tan(beta));
        }
        
    }
    //printf("the value of xmax is: %.12Lf \n",xmax);// same value in mathematica
    //printf("the value of xmin is: %.12Lf \n",xmin);//same value in mathematica
    
    //printf("the value of width is: %.12f \n",width);
    
    
    
    
    for(int n=-20;n<=20;n++)
    {
        for(ex=xmin; ex<xmax; ex+=width/res)
        {
            fc = 2*pi*n*ex/d;
            
            //printf("the value of fc is: %f and n is: %d and ex is %.12f \n",fc,n,ex); //same value in mathematica
            
            ph = -width*thick*chargeratio*pow(e_charge,2)*(2*pi*Coulomb/Plancks)/(vel*(.25*pow(width,2)-pow(ex,2)));
            
            // printf("the value of ph is: %f and n is: %d and ex is %.12f \n",ph,n,ex); //same value in mathematica
            
            
            
            
            ReT[x2pnts((float *)ReT, r, n)][1] += cos(ph+fc);
            ImT[x2pnts((float *)ReT, r, n)][1] += sin(ph+fc);
            
        }
        
        
        
        // for(j=0;j<=41;j++){
        
        //ReT[j][1] += cos(ph+fc);
        //ImT[j][1] += sin(ph+fc);
        
        //printf("the value of ReT is: %f and j is: %d \n",ReT[j][1],j);
        //printf("the value of ImT is: %f and j is: %d \n",ImT[j][1],j);
        //  }
        
        
        
        
    }
    
    
    // testing if the results match with methematica
    for (int i=0; i<=40; i++) {
        //  printf("value of ReT is: %f and j is: %d \n",ReT[i][1],(i-20));////same value in mathematica
        // printf("value of ImT is: %f and j is: %d \n",ImT[i][1],(i-20));////same value in mathematica
        
        
        ReT[i][1] = ReT[i][1]/res;
        ImT[i][1] = ImT[i][1]/res;
        
        //printf("value of ReT is: %0.10f and j is: %d \n",ReT[i][1],(i-20));////same value in mathematica
        //printf("value of ImT is: %0.10f and j is: %d \n",ImT[i][1],(i-20));////same value in mathematica
        
    }
    
    //printf("the value of ReT is: %f and j is: %d \n",ReT[j][1],j);
    
    
    
    
    
    
    
    
    double ix[300][2]={0};
    
    
    for (int i=0; i<300; i++) {
        ix[i][0]= xstart+(i-1)*((xend-xstart)/(xpnts-1));
        //printf("the values of ix[i] and i are: %f\t \t %d\n",ix[i][0],i);//different format, but same values as in the Mathematica code
    }
    
    
    double phix[300][2]={0};
    
    
    for (int i=0; i<300; i++) {
        phix[i][0]= xstart+(i-1)*((xend-xstart)/(xpnts-1));
        //printf("the values of phix[i] and i are: %f\t \t %d\n",phix[i][0],i);//different format, but same values as in the Mathematica code
    }
    
    
    for (int i=0; i<300; i++){
        for (int m1=-lim; m1<=lim; m1++) {
            for (int m2=-lim; m2<=lim; m2++) {
                for (int n1=-lim; n1<=lim; n1++) {
                    for (int n2=-lim; n2<=lim; n2++) {
                        
                        
                        dn =n1-n2;
                        n = ((double)(n1+n2))/2;
                        dm = m1-m2;
                        m = ((double)(m1+m2))/2;
                        
                        
                        a = (x2pnts((float *)ReT, r, m1));
                        b = (x2pnts((float *)ReT, r, m2));
                        c = (x2pnts((float *)ReT, r, n1));
                        d5 = (x2pnts((float *)ReT, r, n2));
                        
                        
                        
                        if (test1==1)
                        {
                            
                            coef = sinc(eta1*pi*m1)+ 0*_Complex_I;
                            coef = coef*(sinc(eta1*pi*m2+ 0*_Complex_I));
                            //printf("the value of coef is: %f \n", coef);
                            
                        }
                        else
                        {
                            
                            
                            coef = ReT[a][1]+ImT[a][1]*_Complex_I;
                            //printf("the value of coef and a are: %0.15f \t %0.15fi \t %d \n", crealf(coef), cimagf(coef), a);
                            //printf("the value of coef and a are:%0.15f \t %0.15fi \t %d \n", cimagf(coef),ImT[a][1], a);
                            
                            coef = coef*((ReT[b][1]-ImT[b][1]*_Complex_I));
                            //printf("the value of coef is: %0.15f \t %0.15fi \t %d \n", ReT[b][1], ImT[b][1], b);
                        }
                        
                        
                        
                        
                        coef = coef*(ReT[c][1] + ImT[c][1]*_Complex_I);
                        
                        //printf("the value of coef is: %f \t %f \n", creal(coef), cimag(coef));
                        
                        coef = coef*(ReT[d5][1] + ImT[d5][1]*_Complex_I);
                        
                        //printf("the value of coef is: %0.10f \t %0.10f \n", creal(coef), cimag(coef));
                        
                        
                        coef=coef*(exp(-pi*pow(((dn*sin(theta)*lambda*(z23))/(d2*el3y)),2)));
                        
                        //printf("the value of coef is: %0.10f \t %0.10f \n", creal(coef), cimag(coef));
                        
                        
                        coef=coef*(exp(-pi*pow((lambda*z23*(dn*cos(theta)+dm*z13/z23)/(d1*el3x)),2)));
                        
                        //printf("the value of coef is: %0.19f \t %0.19fi \n", crealf(coef), cimagf(coef));// results match with mathematica's!
                        
                        
                        
                        
                        
                        
                        
                        
                        if (((__real__ coef)>=cutoff) || ((__imag__ coef)>=cutoff)) {
                            //printf("the value of coef is: %0.10f \t %0.10fi \n", crealf(coef), cimagf(coef));//results match with mathematica's!
                            
                            phi = dn*n*(1-z23/v3x)*pow((cos(theta)),2) + dn*n*(1-z23/v3y)*pow((sin(theta)),2) + dn*m*(1-z13/v3x)*cos(theta);
                            //phi = dn*n*(1-z23/v3x)+ dn*m*(1-z13/v3x)*cos(theta);
                            
                            //printf("the values of phi is: %f \n", phi);
                            
                            
                            phi = phi +(dm*n*(1-z13/v3x)*cos(theta) + dm*m*(z13/z23)*(1-z13/v3x));
                            //printf("the values of phi is: %f \n", phi);
                            
                            phi = phi*(2*pi*lambda*z23/(pow(d1,2)));
                            //printf("the values of phi is: %f \n", phi);
                            
                            phi = phi - (2*pi*dn*G2_x/d2);
                            
                            
                            
                            
                            
                            phix[i][1] = ((phi-(2*pi*(phix[i][0])/d2)*(dn*cos(theta)*(1-z23/v3x) + dm*(1-z13/v3x))));
                            
                            
                            //printf("the values of phix[i][1] is: %f \t %d \n", phix[i][1], i);//values apparently are the same as in mathematica(zeros)
                            
                            
                            
                            
                            ix[i][1] = ix[i][1] + ((((__real__ coef)*cos(phix[i][1]) - (__imag__ coef)*sin(phix[i][1]))*exp(-pi*pow(((phix[i][0]-(lambda*z23/d1)*(n*cos(theta)+m*(z13/z23)))/w3x),2))));
                            
                            
                            
                            
                        }
                        
                    }
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                }
                
            }
            
        }
       //printf("the values of ix[i][1] are:  %0.19f \t %d \n",ix[i][1],i);// the results match with mathematica's!!
      
        
    }
    
    
    double ix3[300]={0};
    double ix4[300]={0};
    
    
    
    for (int i=0; i<300; i++) {
        for (int j=0; j<2; j++) {
            if (j==0) {
                ix3[i]=ix[i][j];
            }
            if (j==1) {
                ix4[i] = ix[i][j];
            }
        }
    }
    
    
    SimplePlot::graph("gp2 graph",ix3,ix4,300);
    
    
  
    
    
}





