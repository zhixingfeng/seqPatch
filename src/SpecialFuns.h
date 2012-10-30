#ifndef SPECIALFUNS_H
#define SPECIALFUNS_H
#define my_pi 3.141592653589793238
#include <math.h>
inline double digamma(double x)
{
    if (x<=0) return sqrt(-1);
    return log(x) - 1/(2*x) - 1/(12*x*x) + 1/(120*pow(x,4)) - 1/(252*pow(x,6));
}

inline double lgamma(double x)
{
    if (x<0.5) return sqrt(-1);
    return (x- 0.5)*log(x) - x + log(2*my_pi)/2+ log( 1+ 1/(12*x) + 1/(288*pow(x,2)) - 139/(51840*pow(x,3)) - 571/(2488320*pow(x,4)) );
}

#endif // SPECIALFUNS_H
