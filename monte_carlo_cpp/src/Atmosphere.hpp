#ifndef ATMOSPHERE_HPP
#define ATMOSPHERE_HPP

#include <string>
#include <vector>

struct AtmosphereLayer
{
    double altitude_m;
    double pressure_pa;
    double temperature_k;
    double molecular_number_density_m3;
    double ozone_number_density_m3;
    double h2o_number_density_m3;
    double no2_number_density_m3;
    double aerosol_extinction_550_m_inv;
    double aerosol_single_scattering_albedo;
    double aerosol_asymmetry;
    double aerosol_scattering_angstrom_exponent;
    double aerosol_absorption_angstrom_exponent;
};

struct AtmosphereState
{
    double altitude_m;
    double pressure_pa;
    double temperature_k;
    double molecular_number_density_m3;
    double ozone_number_density_m3;
    double h2o_number_density_m3;
    double no2_number_density_m3;
    double aerosol_extinction_550_m_inv;
    double aerosol_single_scattering_albedo;
    double aerosol_asymmetry;
    double aerosol_scattering_angstrom_exponent;
    double aerosol_absorption_angstrom_exponent;
};

class Atmosphere
{
public:
    void loadProfileCsv(const std::string &filename);
    void loadLayerData(const std::string &filename);

    const std::vector<AtmosphereLayer> &layers() const { return layers_; }
    bool empty() const { return layers_.empty(); }

    double topOfAtmosphereAltitudeM() const;
    AtmosphereState stateAtAltitude(double altitude_m) const;

    double pressurePa(double altitude_m) const;
    double temperatureK(double altitude_m) const;
    double molecularNumberDensity(double altitude_m) const;
    double ozoneNumberDensity(double altitude_m) const;
    double h2oNumberDensity(double altitude_m) const;
    double no2NumberDensity(double altitude_m) const;
    double aerosolExtinction550(double altitude_m) const;
    double aerosolSingleScatteringAlbedo(double altitude_m) const;
    double aerosolAsymmetry(double altitude_m) const;
    double aerosolScatteringAngstromExponent(double altitude_m) const;
    double aerosolAbsorptionAngstromExponent(double altitude_m) const;

    // Compatibility accessors retained for older code paths.
    double rayleighDensity(double altitude_m) const;
    double aerosolDensity(double altitude_m) const;
    double absorptionCoeff(double altitude_m) const;
    double atmosphericDensity(double altitude_m) const;

private:
    std::vector<AtmosphereLayer> layers_;
};

#endif // ATMOSPHERE_HPP
