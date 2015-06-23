//
//  x2pnts.c
//  
//
//  Created by Arthur Romero on 6/23/15.
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>


int x2pnts(float a[][2], int rows, int columns, int value);

int main(void)
{
    
    double d = 0.0000001;// period of grating
    double e_charge = 0.00000000000000000016021765;
    double e_mass = 0.00000000000000000000000000000091093819;
    double Plancks = 0.0000000000000000000000000000000006626068;
    double Coulomb = 0.00000000898755179;
    double pi = 3.14159265358979;
    double const_e = 2.71828182845905;
    double chargeratio =0;//strength of image charge(units of e (?))
    
    double lambda = 0.00000000001;

    double energy = (1.5*pow(10,-18)/(pow(lambda,2)));
    double width = 4*pow(10,-8); //eta2*d
    double thick = 14*pow(10,-9);
    double wedgeangle = 0;
    double tilt =0;
    
    
    double res = 1000;
    double eta = width/d;
    double vel = pow(2*energy*e_charge/e_mass,1/2);
    double alpha = wedgeangle*pi/180;
    double beta = tilt*pi;// since tilt=0, beta=0.
    
    
    double eta1 = .4;//G1 open fraction
    double eta2 = .4;//G2 open fraction
    
    
    double r0 = -4.04;//initial radius of wavefront curvature
    double el0 = 0.000001;// initial coherence width
    double w0 = 0.00003;// initial beam width
    
    double G1_z = 0.000001;
    double G2_z = 1;
    double G2_x = 0.00000005; //d/2;
    
    double theta = 0;

    
    long double xmin;
    long double xmax;
    

    
    
    
    
    double zstart = -0.1;
    double zend = 2.1;
    double xstart = -200*pow(10,-6);
    double xend = 200*pow(10,-6);
    double ystart = -0.11*pow(10,-3);
    double yend = 0.11*pow(10,-3);
    
    
    
    
    double xpnts = 300;
    double ypnts = 300;
    double zpnts = 300;
    
    
    
    double zres = (zend-zstart)/zpnts;//
    
    
    double w1;
    double r1;
    double el1;
    double el2;
    
    
    float fc;
    float ph;
    float ex;
    int j;

    
    
    
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
        for(ex=xmin; ex<xmax; ex+=width/res)
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
        
        printf("value of ReT is: %f and j is: %d \n",ReT[i][1],(i-20));////same value in mathematica
        //printf("value of ImT is: %f and j is: %d \n",ImT[i][1],(i-20));////same value in mathematica
        
    }

    float index, value;
    int rows,columns;
 
    
    value = -16;/////example
    index = x2pnts(ReT, 41, 2, value);
    int i2 = index;
    if (index == -1)
    {
        printf("The value %f was not found.\n", value);
        
    }
    else
    {
        printf("The value %f was found at %f\n", value, index);
         printf(" the equivalent is: %f\n", ReT[i2-1][1]);
    }
    
    value = 13;//////example
    index = x2pnts(ReT, 41,2, value);
    int i3 = index;
    if (index == -1)
    {
        printf("The value %f was not found.\n", value);
    }
    else
    {
        printf("The value %f was found at %f\n", value, index);
        printf(" the equivalent is: %f\n", ReT[i3-1][1]);
    }
    
    value = 0;////////example
    index = x2pnts(ReT, 41, 2, value);
    int i4 = index;
    if (index == -1)
    {
        printf("The value %f was not found.\n", value);
    }
    else
    {
        printf("The value %f was found at %f\n", value, index);
        printf(" the equivalent is: %f\n", ReT[i4-1][1]);
    }
}

int x2pnts(float a[][2], int rows, int columns, int value)
{
    int i;
    for (i=0; i<rows; i++)
    {
        for (int j=0; j<columns; j++)
        {
            
            
            if (a[i][j] == value)
            {
                return(i+1);  /* it was found */
               
            }
        }
    }
    return(-1);  /* if it was not found */
}

