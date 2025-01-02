// AccelerationTechniques.cpp
#include "AccelerationTechniques.hpp"
#include <random>

// example stubs
static std::mt19937 rngAcc(5678);
static std::uniform_real_distribution<double> uniAcc(0.0,1.0);

double russianRoulette(double &photonWeight, double threshold)
{
    if(photonWeight < threshold)
    {
        double randv = uniAcc(rngAcc);
        if(randv < 0.5){
            // kill
            photonWeight = 0.0;
            return 0.0;
        } else {
            // boost
            photonWeight *= 2.0;
            return photonWeight;
        }
    }
    return photonWeight;
}

void photonSplitting(/*...*/)
{
    // TODO: if region is very scattering, create multiple photons with partial weight
}
