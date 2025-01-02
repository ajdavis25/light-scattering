// PhaseFunctions.cpp
#include "PhaseFunctions.hpp"
#include <cmath>
#include <algorithm>

double rayleighPhase(double cosTheta)
{
    return 0.75 * (1.0 + cosTheta*cosTheta);
}

double henyeyGreenstein(double cosTheta, double g)
{
    double g2 = g*g;
    double denom = 1.0 + g2 - 2.0*g*cosTheta;
    return (1.0 - g2)/std::pow(denom, 1.5);
}

double combinedPhase(double cosTheta, double fracRay, double g)
{
    // Weighted sum of Rayleigh & HG
    double pRay = rayleighPhase(cosTheta);
    double pHG  = henyeyGreenstein(cosTheta, g);
    return fracRay*pRay + (1.0 - fracRay)*pHG;
}
