/******************************************************
 * MonteCarloSky_Advanced.cpp
 *
 * Compile example (Linux or MinGW):
 *    g++ -O3 -std=c++17 MonteCarloSky_Advanced.cpp -o MonteCarloSky
 *
 * Run:
 *    ./MonteCarloSky
 *
 * Key Upgrades over the simple version:
 * 1) Large photon count (1e7) for smoother results
 * 2) Simple absorption function
 * 3) Rayleigh + Henyey-Greenstein aerosol scattering
 * 4) Lambertian Earth surface reflection (albedo)
 * 5) Mixture of Rayleigh & aerosol fraction, altitude-dependent
 *
 * This code is still a demonstration: no spectral 
 * loops, no polarization, no advanced molecular lines, etc.
 ******************************************************/

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>

// If M_PI isn't defined by <cmath>, define our own
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

//--------------------------------------------------
// 1) MODEL CONSTANTS
//--------------------------------------------------
static const double R_EARTH   = 6.371e6;    // Earth radius (m)
static const double ATM_HEIGHT= 1.0e5;      // ~100 km atmosphere
static const double R_TOA     = R_EARTH + ATM_HEIGHT;
static const double H_SCALE   = 8.0e3;      // scale height for Rayleigh
static const double SIGMA_RAY_0 = 1e-5;     // baseline Rayleigh near ground

// "Aerosol fraction" param: we model aerosol with H-G phase function
// fraction for scattering vs. Rayleigh
static const double AEROSOL_FRAC_0 = 0.2;   // fraction at ground (example)
static const double H_SCALE_AER    = 2.0e3; // aerosols might have smaller scale height

// Henyey-Greenstein g-factor for forward-lobe aerosol
static const double G_AEROSOL = 0.7;

// Simple absorption, e.g. from O2 or other gases:
static const double ABS_0 = 2e-5;           // baseline absorption at ground
static const double H_SCALE_ABS = 7.0e3;    // scale height for absorption

// Lambertian surface reflection:
static const double ALBEDO = 0.2;           // 20% reflectivity

// How many photons (big for better stats)
static const long long NUM_PHOTONS = (long long)1e7;

// Binning (zenith up to 90°, azimuth 0..360 in 2° increments => 180 bins)
static const int N_ZEN = 90;
static const int N_AZI = 180;

//--------------------------------------------------
// 2) Random Generator
//--------------------------------------------------
static std::mt19937 rng(12345);
static std::uniform_real_distribution<double> uni(0.0,1.0);
inline double rand01() { return uni(rng); }

//--------------------------------------------------
// 3) Atmosphere Profiles
//--------------------------------------------------
/**
 * Exponential Rayleigh density: 
 *   rho_ray(r) ~ exp(-(r - R_EARTH)/H_SCALE)
 * We'll define a "rayleighCoefficient(r)" that includes SIGMA_RAY_0.
 */
double rayleighCoefficient(double r) {
    double alt = r - R_EARTH;
    if (alt < 0.0) return 0.0;
    return SIGMA_RAY_0 * std::exp(-alt / H_SCALE);
}

/**
 * Exponential aerosol coefficient:
 *   aerosols often have different scale height
 */
double aerosolCoefficient(double r) {
    double alt = r - R_EARTH;
    if (alt < 0.0) return 0.0;
    // some baseline * exp(-alt/H_SCALE_AER)
    // We'll assume "AEROSOL_FRAC_0 * SIGMA_RAY_0" 
    double base = AEROSOL_FRAC_0 * SIGMA_RAY_0; 
    return base * std::exp(-alt / H_SCALE_AER);
}

/**
 * Absorption coefficient (simple exponential).
 *   Could represent O2, O3, etc. 
 */
double absorptionCoefficient(double r) {
    double alt = r - R_EARTH;
    if (alt < 0.0) return 999999.0; // below ground => large => no penetration
    return ABS_0 * std::exp(-alt / H_SCALE_ABS);
}

/**
 * Combined scattering coefficient at radius r:
 *   sigma_scatter = rayleigh(r) + aerosol(r)
 */
double scatteringCoefficient(double r) {
    return rayleighCoefficient(r) + aerosolCoefficient(r);
}

//--------------------------------------------------
// 4) Phase Functions
//--------------------------------------------------
/**
 * Rayleigh phase function (unpolarized single scatter):
 *   p_ray(cosTheta) ~ (3/4)*(1 + cos^2Theta)
 */
inline double phaseRayleigh(double cosTheta) {
    return 0.75 * (1.0 + cosTheta*cosTheta);
}

/**
 * Henyey-Greenstein phase function:
 *   p_HG(cosTheta) = (1 - g^2) / ( (1 + g^2 - 2g cosTheta)^(3/2) )
 * for g in [0..1)
 */
inline double phaseHG(double cosTheta, double g) {
    double g2 = g*g;
    double denom = 1.0 + g2 - 2.0*g*cosTheta;
    return (1.0 - g2) / std::pow(denom, 1.5);
}

/**
 * Combined Phase Function:
 *   fractionRay = rayleigh(r)/[rayleigh(r) + aerosol(r)]
 *   fractionAero= aerosol(r)/[rayleigh(r) + aerosol(r)]
 * Then total p = fractionRay * p_ray(cosTheta) + fractionAero * p_HG(cosTheta,g)
 */
double phaseFunction(double cosTheta, double r) {
    double ray = rayleighCoefficient(r);
    double aer = aerosolCoefficient(r);
    double total = ray + aer;
    if (total <= 0.0) {
        // fallback
        return 0.0;
    }
    double fracRay = ray / total;
    double fracAer = aer / total;
    double pRay = phaseRayleigh(cosTheta);
    double pHG  = phaseHG(cosTheta, G_AEROSOL);
    return fracRay * pRay + fracAer * pHG;
}

// We will do a naive "rejection sampling" for scattering angles again.
// For better performance, you might do separate sampling for Rayleigh & HG, then 
// randomly choose which event occurs. But let's keep it simpler in one step.

// We'll guess an upper bound for p: Rayleigh max~1.5, HG can have a peak >>1. 
// We'll pick something a bit large, say 3.0 or 4.0. You can refine if needed.
static const double PHASE_BOUND = 4.0; 

//--------------------------------------------------
// 5) Simple Vector & Photon
//--------------------------------------------------
struct Vec3 {
    double x, y, z;
};

inline double dot(const Vec3 &a, const Vec3 &b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
inline double norm(const Vec3 &a){ return std::sqrt(dot(a,a)); }

Vec3 normalize(const Vec3 &v){
    double n = norm(v);
    if(n<1e-15) return {0,0,1};
    return {v.x/n, v.y/n, v.z/n};
}

// Photon structure
struct Photon {
    Vec3 pos;
    Vec3 dir; 
    double weight; // track partial absorption, etc.
};

/**
 * Local coords approach to rotate oldDir into new random scatter direction.
 * We'll sample cosTheta in [-1..1], phi in [0..2pi], do naive rejection 
 * with phaseFunction(cosTheta,r).
 */
Vec3 sampleNewDirection(const Vec3 &oldDir, double r) {
    while(true){
        double mu = 2.0*rand01() - 1.0;   // cosTheta
        double phi= 2.0*M_PI*rand01();
        double testVal = PHASE_BOUND * rand01();
        double pVal = phaseFunction(mu, r);
        if(testVal <= pVal) {
            double sinTheta = std::sqrt(std::max(0.0, 1.0 - mu*mu));
            Vec3 w = oldDir;
            // pick 'up' to define a cross
            double wzAbs = std::fabs(w.z);
            Vec3 up = (wzAbs > 0.9)? Vec3{1,0,0} : Vec3{0,0,1};
            // cross => u
            Vec3 u = {
                up.y*w.z - up.z*w.y,
                up.z*w.x - up.x*w.z,
                up.x*w.y - up.y*w.x
            };
            u = normalize(u);
            // v = w x u
            Vec3 v = {
                w.y*u.z - w.z*u.y,
                w.z*u.x - w.x*u.z,
                w.x*u.y - w.y*u.x
            };
            double cp = std::cos(phi), sp = std::sin(phi);
            Vec3 newDir;
            newDir.x = mu*w.x + sinTheta*(cp*u.x + sp*v.x);
            newDir.y = mu*w.y + sinTheta*(cp*u.y + sp*v.y);
            newDir.z = mu*w.z + sinTheta*(cp*u.z + sp*v.z);
            return normalize(newDir);
        }
    }
}

//--------------------------------------------------
// 6) Bins for exit directions
//--------------------------------------------------
static std::vector<std::vector<double>> exitRadiance; 

//--------------------------------------------------
// 7) Earth Reflection
//--------------------------------------------------
/**
 * We do a Lambertian reflection if photon hits r < R_EARTH.
 * That means we reflect with weight *= ALBEDO, random direction upward hemisphere.
 * If we want to keep scattering, we place it at r=R_EARTH (approx).
 */
Vec3 randomUpwardDirection() {
    // sample cosTheta in [0..1], phi in [0..2pi]
    double mu = rand01();  // cosTheta
    double phi= 2.0*M_PI*rand01();
    double sinTheta = std::sqrt(std::max(0.0,1.0 - mu*mu));
    // standard spherical:
    double x = sinTheta*std::cos(phi);
    double y = sinTheta*std::sin(phi);
    double z = mu; 
    // That points "up" if we define +z as radial outward from Earth center
    return {x,y,z};
}

//--------------------------------------------------
// 8) Launch Photons
//--------------------------------------------------
Photon launchPhoton(){
    // Place top center, direction downward
    Photon p;
    p.pos = {0.0, 0.0, R_TOA};
    p.dir = {0.0, 0.0, -1.0};
    p.weight=1.0;
    return p;
}

//--------------------------------------------------
// 9) Recording
//--------------------------------------------------
void recordExit(const Photon &p){
    // direction => (zen,azi)
    double r = norm(p.pos); // should be > R_TOA
    double cosZen = p.dir.z; 
    if(cosZen>1.0)  cosZen=1.0;
    if(cosZen<-1.0) cosZen=-1.0;
    double zenDeg = std::acos(cosZen)*180.0/M_PI;
    double aziDeg = std::atan2(p.dir.y, p.dir.x)*180.0/M_PI;
    if(aziDeg<0.0) aziDeg+=360.0;

    int iZen = (int)std::floor(zenDeg);
    if(iZen<0) iZen=0;
    if(iZen>=N_ZEN) iZen=N_ZEN-1;
    double aBin = aziDeg/2.0;
    int iAzi = (int)std::floor(aBin);
    if(iAzi<0) iAzi=0;
    if(iAzi>=N_AZI) iAzi=N_AZI-1;

    // Accumulate photon weight
    exitRadiance[iZen][iAzi] += p.weight;
}

//--------------------------------------------------
// 10) Tracing photons
//--------------------------------------------------
void tracePhoton(Photon &p){
    while(true){
        double r = norm(p.pos);
        if(r>R_TOA){
            // escapes top
            recordExit(p);
            return;
        }
        if(r< R_EARTH){
            // hits Earth => reflect or absorb 
            // Lambertian reflection with albedo
            p.weight *= ALBEDO; 
            if(p.weight<1e-12){
                // kill if too small
                return;
            }
            // place photon just above surface
            Vec3 dirUp = randomUpwardDirection();
            // scale "pos" to R_EARTH
            double factor = R_EARTH / r; 
            p.pos.x *= factor; 
            p.pos.y *= factor; 
            p.pos.z *= factor; 
            p.dir = dirUp;
            continue;
        }
        // local scattering and absorption
        double sigma_scat = scatteringCoefficient(r);
        double sigma_abs  = absorptionCoefficient(r);
        double sigma_tot  = sigma_scat + sigma_abs;
        if(sigma_tot<1e-15){
            // negligible atmosphere => escapes effectively
            recordExit(p);
            return;
        }

        // free path
        double s = -std::log(rand01())/sigma_tot;
        // move
        p.pos.x += s*p.dir.x;
        p.pos.y += s*p.dir.y;
        p.pos.z += s*p.dir.z;

        // partial absorption
        double absorpProb = sigma_abs/sigma_tot;
        // kill fraction or reduce weight stochastically
        if(rand01() < absorpProb){
            // absorbed => done
            return;
        }
        // scatter => sample new direction
        // new direction depends on local r
        double r2 = norm(p.pos);
        p.dir = sampleNewDirection(p.dir, r2);

        // Optionally track multi-scattering weight if you want
        // p.weight stays the same if we do purely scattering (no weighting).
        // or you can do p.weight *= (sigma_scat / sigma_tot) for "explicit absorption".
        // We'll keep it simpler here.
    }
}

//--------------------------------------------------
// 11) MAIN
//--------------------------------------------------
int main(){
    // Prepare bins
    exitRadiance.resize(N_ZEN);
    for(int iz=0; iz<N_ZEN; iz++){
        exitRadiance[iz].resize(N_AZI, 0.0);
    }

    // Possibly enable OpenMP if you want parallel:
    // #pragma omp parallel for 
    for(long long i=0; i<NUM_PHOTONS; i++){
        Photon p = launchPhoton();
        tracePhoton(p);
    }

    // Print partial results
    std::cout << std::fixed << std::setprecision(6);
    for(int iz=0; iz<N_ZEN; iz+=10){
        for(int ia=0; ia<N_AZI; ia+=20){
            double val = exitRadiance[iz][ia]/double(NUM_PHOTONS);
            double zenMid = iz+0.5;
            double aziMid = (ia*2.0 + 1.0);
            std::cout<<"Zen="<<zenMid<<", Azi="<<aziMid
                     <<", Radiance="<< val <<"\n";
        }
    }

    std::cout << "Monte Carlo done.\n";
    return 0;
}
