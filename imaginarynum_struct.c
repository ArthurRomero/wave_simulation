//
//  imaginarynum_struct.c
//  
//
//  Created by Arthur Romero on 6/30/15.
//
//

#include <stdio.h>
typedef struct complex{
    float real;
    float imag;
}complex;
complex add(complex n1,complex n2);
complex prod(complex n1, complex n2);
int main(){
    complex n1,n2,temp,temp2;
    printf("For 1st complex number \n");
    printf("Enter real and imaginary respectively:\n");
    scanf("%f%f",&n1.real,&n1.imag);
    printf("\nFor 2nd complex number \n");
    printf("Enter real and imaginary respectively:\n");
    scanf("%f%f",&n2.real,&n2.imag);
    temp=add(n1,n2);
    printf("Sum=%.1f+%.1fi\n",temp.real,temp.imag);
    
    temp2 = prod(n1,n2);
    printf("Sum=%.1f+%.1fi\n",temp2.real,temp2.imag);
    return 0;
}
complex add(complex n1,complex n2){
    
    
    complex temp;
    temp.real=n1.real+n2.real;
    temp.imag=n1.imag+n2.imag;
    return(temp);
}

complex prod(complex n1, complex n2){
    
    
    complex temp;
    temp.real = (n1.real*n2.real)-(n1.imag*n2.imag);
    temp.imag=(n1.real*n2.imag)+(n1.imag*n2.real);
    return(temp);
}