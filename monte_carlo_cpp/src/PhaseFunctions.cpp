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

double miePhase(double cosTheta, double wavelength_nm)
{
    // TODO: read from a LUT or do some approximate formula 
    // For demonstration, let's do a placeholder forward-lobe:
    double factor = (550.0 / wavelength_nm);
    double peak = 5.0 * factor;  // bigger forward peak if shorter wavelength?
    double val = 1.0 + peak * cosTheta;
    return std::max(0.0, val);
}

double combinedPhase(double cosTheta, double fracRay, double fracMie)
{
    // ignoring normalization for brevity
    double pRay = rayleighPhase(cosTheta);
    double pMie = miePhase(cosTheta, 550.0/*placeholder*/);
    return fracRay * pRay + fracMie * pMie;
}
