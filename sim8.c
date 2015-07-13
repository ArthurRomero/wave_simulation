//
//  sim8.c
//  
//
//  Created by Arthur Romero on 7/10/15.
//
//

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>






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



void gp0(double z,double r0,double el0, double w0);//prototype

double * gp1(double z,double r0,double el0, double w0);//prototype

void gp2(double z12,double z23, double mytheta, double el1x, double w1x, double r1x, double el1y, double w1y, double r1y, double G2_x);//prototype

int x2pnts(float *arr, int r, int value);//prototype













int main( int argc, char *argv[] )
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
    
    

    
    
    w1=w(G1_z,r0,el0,w0);//gp1
    //printf("the value of w1 is:  %.12f\n",w1);// matches with Mathematica value of w1
    
    r1 = v(G1_z,r0,el0,w0);//gp1
    //printf("the value of r1 is:  %.12f\n",r1);//matches with Mathematica value of r1
    
    
    el1=el(G1_z,r0,el0,w0);//gp1
    //printf("the value of el1 is:  %.12f\n",el1); //matches with Mathematica value of el1
    
    el2 = el(zstart-G1_z+zres*100, r1, el1, w1);//gp2
    // printf("the value of el2 is:  %.12f\n",el2);// matches with mathematica value
    
    
    
    
    
    double *x1;
    x1 = gp1(zstart - G1_z + zres*100,r1,el1,w1);
    
    
    for (int i=0; i<300; i++) {
         printf("*(x + [%d]) : %f\n", i, (*((x1+(i)*2) + 1)) );
    }
    
    printf("test");
    
    

    
    return 0;
    
    
    
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






void gp0 (double z,double r0,double el0, double w0)
{
    double w1;
    double ix[300][2]={0};
    
    
    for (int i=0; i<300; i++) {
        ix[i][0]= xstart+(i-1)*((xend-xstart)/(xpnts-1));
        ix[i][1]=exp(-pi*pow((ix[i][0]/w1),2));
        //printf("the values of ix[i] and i are: %f\t and %d\n",ix[i][1],i);//
        
    }
    
    
}






double * gp1(double z,double r0,double el0, double w0)
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
    
    return ix;
    
    
    
}


