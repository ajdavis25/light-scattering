// Atmosphere.cpp
#include "Atmosphere.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

void Atmosphere::loadLayerData(const std::string &filename)
{
    // TODO: read text or netCDF file 
    // For now, just define some layers manually:
    layers.clear();
    Layer l1; 
    l1.altBottom=0;   
    l1.altTop=10000; 
    l1.densityRay=1.0;       
    l1.densityAero=0.2;      
    l1.absorptionO3=0.1; 
    layers.push_back(l1);

    Layer l2; 
    l2.altBottom=10000; 
    l2.altTop=20000; 
    l2.densityRay=0.5;
    l2.densityAero=0.05; 
    l2.absorptionO3=0.02;
    layers.push_back(l2);

    // etc...
    // If 'filename' is non-empty, parse it. 
    // This is just a placeholder.
}

double Atmosphere::rayleighDensity(double alt) const
{
    // search layers
    for(const auto &L : layers)
    {
        if(alt>=L.altBottom && alt<L.altTop)
        {
            // simple approach:
            return L.densityRay * std::exp(-(alt - L.altBottom)/8000.0);
        }
    }
    return 0.0;
}

double Atmosphere::aerosolDensity(double alt) const
{
    for(const auto &L : layers)
    {
        if(alt>=L.altBottom && alt<L.altTop)
        {
            // or any advanced logic
            return L.densityAero * std::exp(-(alt - L.altBottom)/2000.0);
        }
    }
    return 0.0;
}

double Atmosphere::absorptionCoeff(double alt) const
{
    for(const auto &L : layers)
    {
        if(alt>=L.altBottom && alt<L.altTop)
        {
            // e.g. scaled by altitude
            double scale = std::exp(-(alt - L.altBottom)/7000.0);
            return L.absorptionO3 * scale * 1e-5;
        }
    }
    return 0.0;
}

double Atmosphere::atmosphericDensity(double alt) const
{
    // for demonstration, let's define "atmosphericDensity"
    // as the sum of Rayleigh + aerosol densities 
    // (though physically, these might not strictly "sum" 
    // in a real sense, but this is your choice).
    double ray = rayleighDensity(alt);
    double aer = aerosolDensity(alt);
    return ray + aer;
}
