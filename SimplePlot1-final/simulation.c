//
//  simulation.c
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

int rows = 300;
int rowsT = 41;
double (*q)[2];


//prototype functions:

double zp(double z, double v);// prototype

double w(double z,double r0, double el0, double w0, double energy);// prototype

double el(double z, double r0, double el0, double w0,double energy);// prototype

double v(double z,double r0, double el0, double w0, double energy);// prototype

double sinc(double x);//prototype

double (*ReTgenerator(double ReT[], double energy))[2];//prototype

double (*ImTgenerator(double ReT[], double energy))[2];//prototype

void gp0(double z,double r0,double el0, double w0,double energy);//prototype

void gp1(double z,double r0,double el0, double w0,double energy);//prototype

void gp2(double z12,double z23, double mytheta, double el1x, double w1x, double r1x, double el1y, double w1y, double r1y, double G2_x, double energy); // prototype

int x2pnts(double *arr, int r, int value);//prototype




int main( )
{

    
    double energy = ((1.5*pow(10,-18))/(pow(0.00000000001,2)))*(1);
    
    double zres = (zend-zstart)/zpnts;


    double w1=w(G1_z,r0,el0,w0,energy);
    double r1 = v(G1_z,r0,el0,w0,energy);
    double el1=el(G1_z,r0,el0,w0,energy);
   
    
    gp0( zstart + 1*zres,r0,el0,w0, energy);
    
    gp1(zstart-G1_z+zres*100,r1,el1,w1, energy);
    
    gp2(G2_z-G1_z,zstart+0*zres,theta,el1,w1,r1,el1,w1,r1,G2_x, energy);
   
}






double zp(double z, double v)
{
    
    double zp;
    zp = (v*z)/(z+v);
    return(zp);
    
}



double w(double z,double r0, double el0, double w0, double energy)
{
    double lambda = sqrt((1.5*pow(10,-18))/(energy));
    double w;
    w = (el0)*(fabs((z)/(zp(z,r0)))*((sqrt((1+(pow(((lambda*zp(z,r0))/(el0*w0)),2)))))));
    
    return(w);
    
}

double el(double z,double r0, double el0, double w0, double energy)
{
    double lambda = sqrt((1.5*pow(10,-18))/(energy));
    double w;
    w = (el0)*(fabs((z)/(zp(z,r0)))*((sqrt((1+(pow(((lambda*zp(z,r0))/(el0*w0)),2)))))));
    
    return(w);
    
}

double v(double z,double r0, double el0, double w0, double energy)
{
    
    double lambda = sqrt((1.5*pow(10,-18))/(energy));
    double v;
    v=(z)/(1-zp(z,r0)/(z*(1+pow(((lambda*zp(z,r0)/(el0*w0))),2))));
    
    return(v);
}



double sinc(double x)// this function avoids a division by 0
{
    
    double sinc;
    
    if (x==0)
    {
        sinc = 1;
    }
    else
    {
        sinc = (sin(x))/x;// this function avoids a division by 0
    }
    
    return(sinc);
}



int x2pnts(int value, int *arr)//finds a value in an array and returns its respective index(ix[i][0])
{
    
    
    for (int i = 0; i < rowsT; i++)
    {
        
        if (value == *(arr+i))
        {
            
            return(i);
            
        }
        
    }
    printf("Error! n or m value does not match with any value of x2pnt array\n");
    
    
}

double (*ReTgenerator(double ReT[], double energy))[2]//it defines all the elements of the array ReT
{
    
    double eta = width/d;
    double vel = pow(2*energy*e_charge/e_mass,1/2);
    double alpha = wedgeangle*pi/180;
    double beta = tilt*pi;
    
    double xmin;
    double xmax;
    
    float fc;
    float ph;
    float ex;
    int j;
    
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
    
    for(int n=-((rowsT-1)/2);n<=((rowsT-1)/2);n++)
    {
        for(ex=xmin; ex<xmax; ex+=width/res)
        {
            fc = 2*pi*n*ex/d;
            
            ph = -width*thick*chargeratio*pow(e_charge,2)*(2*pi*Coulomb/Plancks)/(vel*(.25*pow(width,2)-pow(ex,2)));
            
            j=n+((rowsT-1)/2);
            
            ReT[j] += cos(ph+fc);
            
        }
        
        
        
    }
    
    for (int i=0; i<rowsT; i++) {
        
        ReT[i] = ReT[i]/res;
        
    }
    
    
    
}

double (*ImTgenerator(double ImT[], double energy))[2]//it defines all the elements of the array ImT:
{
    double eta = width/d;
    double vel = pow(2*energy*e_charge/e_mass,1/2);
    double alpha = wedgeangle*pi/180;
    double beta = tilt*pi;
    
    double xmin;
    double xmax;
    
    float fc;
    float ph;
    float ex;
    int j;
    
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
    
    
    
    
    
    for(int n=-((rowsT-1)/2);n<=((rowsT-1)/2);n++)//n=((rowsT-1)/2): for example, if rowsT=41, n=20
    {
        for(ex=xmin; ex<xmax; ex+=width/res)
        {
            fc = 2*pi*n*ex/d;
            
            ph = -width*thick*chargeratio*pow(e_charge,2)*(2*pi*Coulomb/Plancks)/(vel*(.25*pow(width,2)-pow(ex,2)));
            
            j=n+((rowsT-1)/2);
            
            ImT[j] += sin(ph+fc);
            
        }
        
        
        
    }
    
    for (int i=0; i<rowsT; i++) {
        
        ImT[i] = ImT[i]/res;
        
    }
    
    
    
}

void gp0(double z,double r0,double el0, double w0, double energy)//it defines the behavior of the wave before the first grating.
{
    
    double jj;
    double w1 = w(z,r0,el0,w0,energy);
    
    double ix[300][2]={0};
    
    
    for(int i=0; i<rows; i++)
    {
        ix[i][0]= xstart+(i)*((xend-xstart)/(xpnts-1));
        
        jj =pow((ix[i][0]/w1),2);// I had to break it appart to make it work
        
        ix[i][1]=exp(-(pi*jj));
     
        printf("the value of ix[i][1] is %f\n",ix[i][1]);
        
    }
    double ix1[300]={0};
    double ix2[300]={0};
    
    
    
    for (int i=0; i<rows; i++) {
        for (int j=0; j<col; j++) {
            if (j==0) {
                ix1[i]=ix[i][j];
            }
            if (j==1) {
                ix2[i] = ix[i][j];
            }
        }
    }
    
    
    SimplePlot::graph("gp1 graph",ix1,ix2,rows);


}


void gp1(double z12,double r1,double el1, double w1, double energy)//it defines the behavior of the wave between the first and second gratings.
{
    double coef;
    double cutoff=pow(10,-3);
    double lim=5;
    
    
    double lambda = sqrt((1.5*pow(10,-18))/(energy));
    

    double w2=w(z12,r1,el1,w1,energy);
    double r2 = v(z12,r1,el1,w1,energy);
    double el2 = el(z12, r1, el1, w1,energy);
    
    int pos[41]={0};
    for (int i=0; i<rowsT; i++)
    {
        pos[i]=i-((rowsT-1)/2);
    }
    
    
    
    double ReT[41]{0};
    
    q = ReTgenerator(ReT,energy);
    
    double ImT[41]={0};
    
    q = ImTgenerator(ImT,energy);
    
    
    double ix[300][2]={0};
    for (int i=0; i<rows; i++) {
        ix[i][0]= xstart+(i)*((xend-xstart)/(xpnts-1));
        
    }
    
    
    
    for (int i=0; i<rows; i++) {
        for (int n=-lim; n<=lim; n++)
        {
            
            for (int m=-lim; m<=lim; m++)
            {
                double dn =n-m;
                double dm = (m+n)/2;
                double test1=0;
                
                
                
                
                
                
                
                
                
                if (test1==1)
                {
                    coef = sinc(eta1*pi*n)*(sinc(eta1*pi*m)*pow((eta1), 2));
                    
                }
                else
                {
                    coef =  ReT[x2pnts(n, (int *)pos)]*ReT[x2pnts(m,(int *)pos)]+ImT[x2pnts(n,(int *)pos)]*ImT[x2pnts(m,(int *)pos)];//
                    
                    
                }
                
                
                coef = coef*exp(-pi*(dn*(lambda*z12/(pow(d*el2,2)))));// added isfinite macro in order to avoid inf values
                
                if (isfinite(coef)==0)
                {
                    coef=0;
                }
                
                
                if (coef>=cutoff)
                {
                    
                    ix[i][1] = ix[i][1] +  (coef*exp(-pi*pow(((ix[i][0]-dm*lambda*z12/d)/w2),2)*cos(2*pi*(dn/d)*(ix[i][0]-dm*lambda*z12/d)*(1-z12/r2))));
                    
                    
                    
                    continue;
                    
                }
            }
            
            
            
            
            
            
            
            
            
        }
        printf("the values of ix[i][1] are:  %0.19f \t %d \n",ix[i][1],i);
        
        
    }

    
    
    double ix1[300]={0};
    double ix2[300]={0};
    
    
    
    for (int i=0; i<rows; i++) {
        for (int j=0; j<col; j++) {
            if (j==0) {
                ix1[i]=ix[i][j];
            }
            if (j==1) {
                ix2[i] = ix[i][j];
            }
        }
    }
    
    
    SimplePlot::graph("gp1 graph",ix1,ix2,rows);
    
    
    
    
}

void gp2(double z12,double z23, double mytheta, double el1x, double w1x, double r1x, double el1y, double w1y, double r1y, double G2_x, double energy)//defines the behavior of the wave after the second gratings
{
    
    double lambda = sqrt((1.5*pow(10,-18))/(energy));

    double theta = pi*mytheta/180;
    double d1=d;
    double d2=d;
    double z13 = z12+z23;
    
    
    double phi =0;
    
    double cutoff = 0.001;
    
    double lim =5;
    
    double _Complex coef;
    
    double dn = 0;
    double dm =0;
    double m=0;
    double n=0;
    int a5 =0;
    int b  =0;
    int c5 =0;
    int d5=0;
    double test1=0;
    
    
    double el3x = el(z13, r1x, el1x, w1x, energy);//G2z - G1z + zstart + 0*zres, r1, el1, w1
    double w3x = w(z13,r1x,el1x,w1x, energy);
    double v3x = v(z13,r1x,el1x,w1x, energy);
    double el3y = el(z13,r1y,el1y, w1y, energy);
    double w3y = w(z13,r1y,el1y,w1y, energy);
    double v3y = v(z13,r1y,el1y,w1y, energy);
    
    
    
    int pos[41]={0};
    for (int i=0; i<rowsT; i++)
    {
        pos[i]=i-((rowsT-1)/2);
    }
    
    
    
    
    
    double ReT[41]={0};
    
    q = ReTgenerator(ReT,energy);
    
    double ImT[41]={0};
    
    q = ImTgenerator(ImT,energy);
    
    
    
    double ix[300][2]={0};
    
    for (int i=0; i<rows; i++) {
        ix[i][0]= xstart+(i)*((xend-xstart)/(xpnts-1));
        
    }
    
    
    double phix[300][2]={0};
    
    
    for (int i=0; i<rows; i++) {
        phix[i][0]= xstart+(i)*((xend-xstart)/(xpnts-1));
        
    }
    
    
    for (int i=0; i<rows; i++){
        for (int m1=-lim; m1<=lim; m1++) {
            for (int m2=-lim; m2<=lim; m2++) {
                for (int n1=-lim; n1<=lim; n1++) {
                    for (int n2=-lim; n2<=lim; n2++) {
                        
                        
                        dn =n1-n2;
                        n = ((double)(n1+n2))/2;
                        dm = m1-m2;
                        m = ((double)(m1+m2))/2;
                        
                        
                        a5 = (x2pnts(m1, (int *)pos));
                        b  = (x2pnts(m2, (int *)pos));
                        c5 = (x2pnts(n1, (int *)pos));
                        d5 = (x2pnts(n2, (int *)pos));
                        
                        
                        
                        if (test1==1)
                        {
                            
                            coef = sinc(eta1*pi*m1)+ 0*_Complex_I;
                            coef = coef*(sinc(eta1*pi*m2+ 0*_Complex_I));
                            
                            
                        }
                        else
                        {
                            
                            
                            coef = ReT[a5]+ImT[a5]*_Complex_I;
                            
                            coef = coef*((ReT[b]-ImT[b]*_Complex_I));
                            
                        }
                        
                        
                        
                        
                        coef = coef*(ReT[c5] + ImT[c5]*_Complex_I);
                        
                        coef = coef*(ReT[d5] + ImT[d5]*_Complex_I);
                        
                        coef=coef*(exp(-pi*pow(((dn*sin(theta)*lambda*(z23))/(d2*el3y)),2)));
                        
                        coef=coef*(exp(-pi*pow((lambda*z23*(dn*cos(theta)+dm*z13/z23)/(d1*el3x)),2)));
                        
                        
                        
                        
                        
                        
                        if (((__real__ coef)>=cutoff) || ((__imag__ coef)>=cutoff)) {
                            
                            phi = dn*n*(1-z23/v3x)*pow((cos(theta)),2) + dn*n*(1-z23/v3y)*pow((sin(theta)),2) + dn*m*(1-z13/v3x)*cos(theta);
                            
                            phi = phi +(dm*n*(1-z13/v3x)*cos(theta) + dm*m*(z13/z23)*(1-z13/v3x));
                            
                            phi = phi*(2*pi*lambda*z23/(pow(d1,2)));
                            
                            phi = phi - (2*pi*dn*G2_x/d2);
                            
                            
                            
                            
                            
                            phix[i][1] = ((phi-(2*pi*(phix[i][0])/d2)*(dn*cos(theta)*(1-z23/v3x) + dm*(1-z13/v3x))));
                            
                            
                            ix[i][1] = ix[i][1] + ((((__real__ coef)*cos(phix[i][1]) - (__imag__ coef)*sin(phix[i][1]))*exp(-pi*pow(((phix[i][0]-(lambda*z23/d1)*(n*cos(theta)+m*(z13/z23)))/w3x),2))));
                            
                            
                            
                            
                        }
                        
                    }
                    
                    
                }
                
            }
            
        }
        printf("the values of ix[i][1] are:  %0.19f \t %d \n",ix[i][1],i);
    }
    
    
    double ix3[300]={0};
    double ix4[300]={0};
    
    
    
    for (int i=0; i<rows; i++) {
        for (int j=0; j<col; j++) {
            if (j==0) {
                ix3[i]=ix[i][j];
            }
            if (j==1) {
                ix4[i] = ix[i][j];
            }
        }
    }
    
    
    SimplePlot::graph("gp2 graph",ix3,ix4,rows);
    
    
  
    
    
}





