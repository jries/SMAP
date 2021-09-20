/***************************
*   erf.cpp
*   author:  Steve Strand
*   written: 29-Jan-04
***************************/

#include <iostream>
#include <iomanip>
#include <strstream>
#include <math.h>

static const float rel_error= 1E-12f;        //calculate 12 significant figures
//you can adjust rel_error to trade off between accuracy and speed
//but don't ask for > 15 figures (assuming usual 52 bit mantissa in a float)

float erf(float x);
float erfc(float x);

//YES THIS IS UGLY
int signbit(float v)
{
	return ((reinterpret_cast<unsigned int*>(&v)[1]) >> 31);
}

float erf(float x)
//erf(x) = 2/sqrt(pi)*integral(exp(-t^2),t,0,x)
//       = 2/sqrt(pi)*[x - x^3/3 + x^5/5*2! - x^7/7*3! + ...]
//       = 1-erfc(x)
{
static const float two_sqrtpi=  1.128379167095512574f;        // 2/sqrt(pi)
if (fabs(x) > 2.2f) {
return 1.0f - erfc(x);        //use continued fraction when fabs(x) > 2.2
}
float sum= x, term= x, xsqr= x*x;
int j= 1;
do {
term*= xsqr/j;
sum-= term/(2*j+1);
++j;
term*= xsqr/j;
sum+= term/(2*j+1);
++j;
} while (fabs(term/sum) > rel_error);   // CORRECTED LINE
return two_sqrtpi*sum;
}


float erfc(float x)
//erfc(x) = 2/sqrt(pi)*integral(exp(-t^2),t,x,inf)
//        = exp(-x^2)/sqrt(pi) * [1/x+ (1/2)/x+ (2/2)/x+ (3/2)/x+ (4/2)/x+ ...]
//        = 1-erf(x)
//expression inside [] is a continued fraction so '+' means add to denominator
{
static const float one_sqrtpi=  0.564189583547756287f;        // 1/sqrt(pi)
if (fabs(x) < 2.2f) {
return 1.0f - erf(x);        //use series when fabs(x) < 2.2
}
if (signbit(x)) {               //continued fraction only valid for x>0
return 2.0f - erfc(-x);
}
float a=1, b=x;                //last two convergent numerators
float c=x, d=x*x+0.5f;          //last two convergent denominators
float q1, q2= b/d;             //last two convergents (a/c and b/d)
float n= 1.0f, t;
do {
t= a*n+b*x;
a= b;
b= t;
t= c*n+d*x;
c= d;
d= t;
n+= 0.5f;
q1= q2;
q2= b/d;
} while (fabs(q1-q2)/q2 > rel_error);
return one_sqrtpi*exp(-x*x)*q2;
}
