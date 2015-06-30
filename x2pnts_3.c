//
//  x2pnts_3.c
//  
//
//  Created by Arthur Romero on 6/26/15.
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

double x2pnts(double *arr, int value);





int main()
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
    
    double zres = (zend-zstart)/zpnts;
    
    double w1;
    double r1;
    double el1;
    double el2;
    
    
    float fc;
    float ph;
    float ex;
    int j;
    
    double ReT[41][2]={{0}};
    for (int i=0; i<=41; i++) {
        ReT[i][0]=i-20;
        
    }

    double ImT[41][2]={{0}};
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
    // printf("the value of xmax is: %.12Lf \n",xmax);// same value in mathematica
    // printf("the value of xmin is: %.12Lf \n",xmin);//same value in mathematica
    //printf("the value of p is: %.12Lf \n",p);
    //printf("the value of width is: %.12f \n",width);
    
    
    
    
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
        
        printf("value of ReT is: %f and j is: %d \n",ReT[i][1],(i-20));////same value in mathematica
        //printf("value of ImT is: %f and j is: %d \n",ImT[i][1],(i-20));////same value in mathematica
        
    }

    int n = -12.0;// examples of n values//
    double A;
    A = x2pnts((double *)ReT,n);
    
    
    printf("the value of n is: %d\n",n);
    printf("the value of A is: %f\n", A);
    
    
    
    return 0;
}


double x2pnts(double *arr, int value)
{
    int i, j;
    int r = 41, c = 2;
    int arr2[41] = {0};
    
    double ans;
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            //printf("%f \n", *((arr+i*n) + j));
            arr2[i] = *((arr+i*c) + j);
            //printf("the values of arr2 are: %d \t %d \n", arr2[i],i);
            
            if (arr2[i] == value) {
                //printf("the value of i is: %d\n",i+1);
                //printf("the equivalent for %d in ReT is the %d th element: %f\n",value,i+1, *((arr+(i)*2) + 1));
                ans = (*((arr+(i)*2) + 1));
                //printf("the value of ans is: %f\n", ans);
                return(ans);
            }
        }
    }
    
    
    
}

