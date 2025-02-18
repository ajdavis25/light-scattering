/*******************************************************
 * MonteCarloDriver.cpp
 *
 * Multi-scattering driver:
 *   - random twilight launch (theta_sun ~ 85±1°, phi_sun ~ 0±5°)
 *   - partial reflection at surface
 *   - real absorption from atmosphere
 *   - finer binning: 0.5° in zen, 1° in azim
 *   - prints only bins with val>0
 *******************************************************/
#include "MonteCarloDriver.hpp"
#include "Atmosphere.hpp"
#include "WavelengthHandling.hpp"
#include "PhaseFunctions.hpp"
#include "SurfaceReflection.hpp"
#include "Polarization.hpp"

#include <random>
#include <cmath>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Earth's geometry
static const double R_EARTH = 6.371e6;
static const double R_TOA   = R_EARTH + 1.0e5;

// Bins for final sky map
static const int N_ZEN=180;   // 0..180 in 0.5° increments
static const int N_AZI=360;   // 0..360 in 1° increments

static std::vector<std::vector<std::vector<double>>> exitRadiance; 

static Atmosphere atm;
static WavelengthManager waveMgr;

// Thread-local RNG
static thread_local std::mt19937 threadRng(12345);
static thread_local std::uniform_real_distribution<double> threadUni(0.0,1.0);

/**
 * launchPhoton: random near-horizon
 *  - we randomize around ~85° ±1°, phi ~0° ±5° to avoid single bin
 */
Photon launchPhoton(double wave_nm, double fluxFrac)
{
    Photon p;
    // position top-of-atmosphere
    p.pos = {0.0, 0.0, R_TOA};

    // define a small random cone around 85°
    double baseThetaDeg = 85.0;
    double coneThetaDeg = 1.0;  // +/- 1° random
    double basePhiDeg   = 0.0;
    double conePhiDeg   = 5.0;  // +/- 5°

    double randTheta = baseThetaDeg + coneThetaDeg * (2.0*threadUni(threadRng)-1.0);
    double randPhi   = basePhiDeg   + conePhiDeg   * (2.0*threadUni(threadRng)-1.0);

    // clamp angles if needed
    if(randTheta<0.0)   randTheta=0.0;
    if(randTheta>180.0) randTheta=180.0;
    // phi can wrap, that’s fine

    double tsun = randTheta*M_PI/180.0;
    double psun = randPhi  *M_PI/180.0;

    p.dir.x = std::sin(tsun)*std::cos(psun);
    p.dir.y = std::sin(tsun)*std::sin(psun);
    p.dir.z = std::cos(tsun);

    // weight
    p.weight       = fluxFrac;
    p.wavelength_nm= wave_nm;
    p.stokes       = initUnpolarized(p.weight);
    return p;
}

/**
 * recordExitLocal: if r>R_TOA => record final direction
 *   Now we do finer bin:
 *   - iZen = floor(2.0*zenDeg), so 0..180 => 0..360
 *   - iAzi = floor(aziDeg), so 0..360 => 0..359
 */
static void recordExitLocal(const Photon &p, int bandIndex,
                            std::vector<std::vector<double>> &localRadiance)
{
    double r = std::sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);
    if(r < R_TOA) return;

    double cosZen = p.dir.z;
    if(cosZen> 1.0) cosZen= 1.0;
    if(cosZen<-1.0) cosZen=-1.0;
    double zenDeg = std::acos(cosZen)*180.0/M_PI;
    double aziDeg = std::atan2(p.dir.y, p.dir.x)*180.0/M_PI;
    if(aziDeg<0.0) aziDeg += 360.0;

    int iZen = (int)std::floor(2.0*zenDeg);  // 0.5° increments
    if(iZen<0)       iZen=0; 
    if(iZen>=N_ZEN)  iZen=N_ZEN-1;

    int iAzi = (int)std::floor(aziDeg);      // 1° increments
    if(iAzi<0)       iAzi=0; 
    if(iAzi>=N_AZI)  iAzi=N_AZI-1;

    localRadiance[iZen][iAzi] += p.weight;
}

/**
 * tracePhoton: multi-scatter loop
 *  - partial reflection => reflectLambertian(0.2)
 *  - absorption => from atm
 */
void tracePhoton(Photon &p, int bandIndex,
                 std::vector<std::vector<double>> &localRadiance)
{
    while(true)
    {
        double r = std::sqrt(p.pos.x*p.pos.x + p.pos.y*p.pos.y + p.pos.z*p.pos.z);

        // if out top
        if(r>R_TOA){
            recordExitLocal(p, bandIndex, localRadiance);
            return;
        }
        // if below Earth => partial reflection
        if(r<R_EARTH){
            auto refl = reflectLambertian(0.2); // 20% albedo
            p.dir = refl.dir;
            p.weight *= refl.weightMultiplier;
            if(p.weight < 1e-12) return;
            double scale = R_EARTH/r;
            p.pos.x *= scale;
            p.pos.y *= scale;
            p.pos.z *= scale;
            continue;
        }

        double alt = r - R_EARTH;
        double sigmaAbs  = atm.absorptionCoeff(alt);
        double sigmaScat = scatteringCoefficient(atm, alt, p.wavelength_nm);
        double sigmaTot  = sigmaScat + sigmaAbs;
        if(sigmaTot<1e-15){
            recordExitLocal(p, bandIndex, localRadiance);
            return;
        }

        // sample free path
        double s = -std::log(threadUni(threadRng))/sigmaTot;
        p.pos.x += s*p.dir.x;
        p.pos.y += s*p.dir.y;
        p.pos.z += s*p.dir.z;

        // partial absorption
        if(threadUni(threadRng) < (sigmaAbs/sigmaTot)){
            return; // absorbed
        }

        // scattering angle
        double dRay  = atm.rayleighDensity(alt);
        double dAero = atm.aerosolDensity(alt);
        double sumD  = dRay + dAero + 1e-20;
        double fracRay = dRay/sumD;

        // sample from combined phase
        while(true){
            double mu      = 2.0*threadUni(threadRng)-1.0;
            double testVal = 4.0*threadUni(threadRng);

            double val = fracRay*rayleighPhase(mu)
                       + (1.0-fracRay)*henyeyGreenstein(mu,0.7);
            if(testVal<=val){
                // optional polarization
                double phi= 2.0*M_PI*threadUni(threadRng);
                double sinT= std::sqrt(std::max(0.0,1.0 - mu*mu));
                p.dir.z= mu;
                p.dir.x= sinT*std::cos(phi);
                p.dir.y= sinT*std::sin(phi);
                break;
            }
        }

        // kill if too low
        if(p.weight<1e-12) return;
    }
}

/**
 * runMonteCarloSimulation:
 *  - large photon count for stats
 *  - random twilight angle in launchPhoton
 *  - partial reflection, real absorption
 *  - half-degree bins in zen, 1° bins in az
 */
void runMonteCarloSimulation()
{
    atm.loadLayerData(""); // or real data
    const auto &bands = waveMgr.getBands();

    // N_ZEN=180 => 0..179
    // N_AZI=360 => 0..359
    exitRadiance.resize(bands.size());
    for(int b=0; b<(int)bands.size(); b++){
        exitRadiance[b].resize(N_ZEN, std::vector<double>(N_AZI,0.0));
    }

    // e.g. 2 million photons per band
    const long long NPHOT= 2000000;

    #pragma omp parallel for
    for(int b=0; b<(int)bands.size(); b++){
        std::vector<std::vector<double>> localRadiance(N_ZEN, std::vector<double>(N_AZI,0.0));

        for(long long i=0; i<NPHOT; i++){
            Photon p = launchPhoton(bands[b].wavelength_nm, bands[b].fluxFraction);
            tracePhoton(p, b, localRadiance);
        }

        #pragma omp critical
        {
            for(int iz=0; iz<N_ZEN; iz++){
                for(int ia=0; ia<N_AZI; ia++){
                    exitRadiance[b][iz][ia] += localRadiance[iz][ia];
                }
            }
        }
    }

    // print final results => only bins with val>0
    for(int b=0; b<(int)bands.size(); b++){
        std::cout<<"\nBand #" << b
                 <<" wave="<<bands[b].wavelength_nm<<" nm\n";
        for(int iz=0; iz<N_ZEN; iz++){
            for(int ia=0; ia<N_AZI; ia++){
                double val = exitRadiance[b][iz][ia];
                if(val>0.0){
                    // convert iz back to actual zenDeg
                    // since iZen=2.0*zenDeg => zenDeg= iZen/2.0
                    double zenDeg= iz*0.5;
                    // iAzi= floor(aziDeg)
                    double aziDeg= (double)ia;

                    std::cout<<" zen="<<zenDeg
                             <<" azi="<<aziDeg
                             <<" => "<< val <<"\n";
                }
            }
        }
    }

    std::cout<<"Monte Carlo simulation done.\n";
}
