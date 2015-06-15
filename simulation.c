//
//  test1.c
//  
//
//  Created by Arthur Romero on 6/5/15.
//
//

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "simulation.h"





double zp(double z, double v);// prototype

double w(double z,double r0, double el0, double w0);// prototype

double el(double z, double r0, double el0, double w0);// prototype

double v(double z,double r0, double el0, double w0);// prototype



double gp0(double z,double r0,double el0, double w0,double ix_);//prototype


int main( void ){
    
    
double e_charge = 0.00000000000000000016021765;
double e_mass = 0.00000000000000000000000000000091093819;
double Plancks = 0.0000000000000000000000000000000006626068;
double Coulomb = 0.00000000898755179;
double pi = 3.14159265;
double const_e = 2.71828183;



double chargeratio =0;//strength of image charge(units of e (?))


    
double lam = pow((150)/(4000),1/2)*pow(10,-10);// 1.936E-10//

double lambda = 0.000000000001;



double eta1 = .5;//G1 open fraction
double eta2 = .5;//G2 open fraction
double d = 0.0000001;// period of grating

double r0 = -4.04;//initial radius of wavefront curvature
double el0 = 0.000001;// initial coherence width
double w0 = 0.00003;// initial beam width

double G1_z = 0.000001;
double G2_z = 1;
double G2_x = d/2;

double theta = 0;


///////////Fourier Components Variables//////////


double energy = (1.5*pow(10,-18)/(pow(lambda,2)));
double width = 4*pow(10,-8); //eta2*d
double thick = 14*pow(10,-9);
double wedgeangle = 0;
double tilt =0;
//double ReT[41][2]={0};
    
    float ReT[41][2]={{0}};
    for (int i=0; i<=41; i++) {
        ReT[i][0]=i-20;
        
    }

    
    
    float ImT[41][2]={{0}};
    for (int i=0; i<=41; i++) {
        ImT[i][0]=i-20;
        
    }
    
    
    
    
double res = 1000;
double eta = width/d;
double vel = pow(2*energy*e_charge/e_mass,1/2);
double alpha = wedgeangle*pi/180;
double beta = tilt*pi;// since tilt=0, beta=0.


    long double xmin;
    long double xmax;
    //long double p;
    
    


    if (beta>=0){
        
        
        
        xmin= width*(1/res - cos(beta)/2);
        //p= width*(1/res - cos(beta)/2);
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
    // printf("the value of xmax is: %.12Lf \n",xmax);// same value in mathematica
   // printf("the value of xmin is: %.12Lf \n",xmin);//same value in mathematica
    //printf("the value of p is: %.12Lf \n",p);
    //printf("the value of width is: %.12f \n",width);
    


        
        
        
        
    float fc;
    float ph;
    float ex;
    int j;
    

        for(int n=-20;n<=20;n++)
        {
            for(ex=xmin; ex<xmax; ex+=width/res)
            {
                fc = 2*pi*n*ex/d;
                
                //printf("the value of fc is: %f and n is: %d and ex is %.12f \n",fc,n,ex); //same value in mathematica

                ph = -width*thick*chargeratio*pow(e_charge,2)*(2*pi*Coulomb/Plancks)/(vel*(.25*pow(width,2)-pow(ex,2)));
                
                //printf("the value of ph is: %f and n is: %d and ex is %.12f \n",ph,n,ex); //same value in mathematica
                
                
                j=n+20;
                
                ReT[j][1] += cos(ph+fc);
                ImT[j][1] += sin(ph+fc);
                
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
        
        //printf("value of ReT is: %f and j is: %d \n",ReT[i][1],(i-20));////same value in mathematica
        //printf("value of ImT is: %f and j is: %d \n",ImT[i][1],(i-20));////same value in mathematica
        
    }
    
   // printf("the value of ReT is: %f and j is: %d \n",ReT[j][1],j);
    
    
    
    
    
    
    
    
    
    
    
 /////////////////////////////////////----Set spatial limits of simulation and initialize waves containg simulation outputs-----//////////////////////////////////////////////
    
    
    
    double zstart = -0.1;
    double zend = 2.1;
    double xstart = -200*pow(10,-6);
    double xend = 200*pow(10,-6);
    double ystart = -0.11*pow(10,-3);
    double yend = 0.11*pow(10,-3);
    
    
    
    
    double xpnts = 300;
    double ypnts = 300;
    double zpnts = 300;
    
    
    
    double zres = (zend-zstart)/zpnts;// matches with Mathematica code value
    
    //printf("the value of zres is: %f\n",zres);
    
    

    
    
    float fmap[150][300];    //
    
    for (int i=0; i<=150; i++) {
        for (int j=0; j<=300; j++) {
            
            fmap[i][j]=0;
   
    //printf("the value of fmap1 is: %f , i = %d \t , j = %d\n",fmap1[i][j],i,j); /// fmap is an array of 150 rows and 300 columns
    }
    }
    
    
    
    float izx[300][300];
    
    for (int i=0; i<=300; i++) {
        for (int j=0; j<=300; j++) {
            
            izx[i][j]=0;
            
            //printf("the value of izx is: %f , i = %d \t , j = %d\n",izx[i][j],i,j); // izx is an array of 300 rows and 300 columns
        }
    }
    

    
    
    
    
    double ix[300][1]={0};
    
    
    for (int i=0; i<=300; i++) {
        ix[i][0]= xstart+(i-1)*((xend-xstart)/(xpnts-1));
       printf("the values of ix[i] and i are: %f\t and %d\n",ix[i][0],i);//different format, but same values as in the Mathematica code
    }
    
    
    double g1[600][1]={0}; // instead of 0, it should be NaN // I do not know the difference
    
    
    for (int i=0; i<=600; i++) {
        g1[i][0]= xstart+(i-1)*((xend-xstart)/(2*xpnts-1));
        //printf("the valueof g1[i][0] at i, %f\t is: %d\n",g1[i][0],i); //different format, but same values as in the Mathematica code
    }
    

    
    
    double g2[600][1]={0}; // instead of 0, it should be NaN // I do not know the difference
    
    
    for (int i=0; i<=600; i++) {
        g2[i][0]= xstart+(i-1)*((xend-xstart)/(2*xpnts-1));
        //printf("the valueof g2[i][0] at i, %f\t is: %d\n",g2[i][0],i);//different format, but same values as in the Mathematica code
    }
 
    
    
    double vis[300][1]={0}; //
    
    
    for (int i=0; i<=300; i++) {
        vis[i][0]= (zstart-G1_z)+(i-1)*((zend-zstart)/(zpnts-1));
        //printf("the valueof vis[i][0] at i, %f\t is: %d\n",vis[i][0],i);//different format, but same values as in the Mathematica code
    }

    
    
    
    
    //lacks ixy
    
    ////////////////////////////////////////////////////////////////
    
    
    ///////Calculating GSM parameter at first grating://///
    
    
    
    
    
    
    //gp0:
    
    
    double w1;
    double r1;
    double el1;
    double ix_;

    
    w1=w(G1_z,r0,el0,w0);
    
    
   //printf("the value of w1 is:  %.12f\n",w1);// matches with Mathematica value of w1
    
 
    
    r1 = v(G1_z,r0,el0,w0);
    
    
   // printf("the value of r1 is:  %.12f\n",r1);//matches with Mathematica value of r1
    
    
    el1=el(G1_z,r0,el0,w0);
    
    
   // printf("the value of el1 is:  %.12f\n",el1); //matches with Mathematica value of el1
    
    //// Range of plotting from Mathematica: x:{-200*10^-6, -170*10^-6}, y:{-10^-6, 10^-6}}
    
    
    //gp0:
    
    
    for (int i=0; i<=300; i++) {
        
        
        ix_= ix[i][1];
        
    ix[i][1]=gp0(zstart + zres*1, r0, el0, w0,ix[i][1]);
        printf("value of ix: %.12f\t and i: %d \n",ix[i][0],i);
    }
 
    
    
    
    
    
    
    
    
   //gp0 = gp0(zstart+zres*1,r0,el0,w0);

    
    
    //gp1:
    
    
    double el2;
    
    
    
    el2 = el(zstart-G1_z+zres*100, r1, el1, w1);
    
    //printf("the value of el2 is:  %.12f\n",el2);// might be wrong
    
    
    
    
    
    
    
    
    
    
    
        
}



double zp(double z, double v)
{

    double zp;
    zp = v*z/(z+v);
    return(zp);
    
}



double w(double z,double r0, double el0, double w0)
{

double w;
double lambda = 0.000000000001;

w = el0*fabs((z)/(((zp(z,r0)))))*pow(1+pow(lambda,2)*pow(zp(z,r0),2)/((pow((el0*w0),2))),1/2);
                    
return(w);

}

double el(double z,double r0, double el0, double w0)
{
    
    double w;
    double lambda = 0.000000000001;
    
    w = el0*fabs((z)/(((zp(z,r0)))))*pow(1+pow(lambda,2)*pow(zp(z,r0),2)/((pow((el0*w0),2))),1/2);
    
    return(w);
    
}

double v(double z,double r0, double el0, double w0)
{
    
    double v;
    
    double lambda = 0.000000000001;
    
    v=(z)/(1-zp(z,r0)/(z*(1+pow(((lambda*zp(z,r0)/(el0*w0))),2))));
    
    return(v);
}



double gp0(double z,double r0,double el0, double w0,double ix_)
{
  double w1;
  double pi = 3.14159265;
    
    w1 = w(z,r0,el0,w0);
    
  
        
        
        
        
    ix_=exp(-pi*pow((ix_/w1),2));
        
   

    return(ix_);

}
