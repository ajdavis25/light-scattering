// Atmosphere.hpp
#ifndef ATMOSPHERE_HPP
#define ATMOSPHERE_HPP

#include <vector>
#include <string>

struct Layer
{
    double altBottom;     // bottom altitude (m)
    double altTop;        // top altitude (m)
    double densityRay;    // Rayleigh reference
    double densityAero;   // aerosol reference
    double absorptionO3;  // or combined absorption
    // possibly temperature, pressure, etc.
};

class Atmosphere
{
public:
    /**
     * Load layering from a file or set default exponentials
     */
    void loadLayerData(const std::string &filename);

    /**
     * Returns Rayleigh density at altitude
     */
    double rayleighDensity(double alt) const;

    /**
     * Returns aerosol density at altitude
     */
    double aerosolDensity(double alt) const;

    /**
     * Returns absorption coefficient at altitude
     */
    double absorptionCoeff(double alt) const;

    /**
     * Returns *total* atmospheric density at altitude,
     * or some other combined measure. You can define
     * how you want this to sum or combine.
     */
    double atmosphericDensity(double alt) const;

private:
    std::vector<Layer> layers;
};

#endif // ATMOSPHERE_HPP
