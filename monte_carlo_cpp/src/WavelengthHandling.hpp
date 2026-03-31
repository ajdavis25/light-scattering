#ifndef WAVELENGTHHANDLING_HPP
#define WAVELENGTHHANDLING_HPP

#include <string>
#include <utility>
#include <vector>

class Atmosphere;

struct SpectralBand
{
    double wavelength_nm;
    double solar_irradiance_w_m2_nm;
    double instrument_response;
    double normalized_weight;
    double rayleigh_cross_section_m2;
    double ozone_cross_section_m2;
    double o2_cross_section_m2;
    double o4_cross_section_m5;
    double h2o_cross_section_m2;
    double no2_cross_section_m2;
};

struct LocalOpticalProperties
{
    double rayleigh_scattering_m_inv;
    double aerosol_scattering_m_inv;
    double aerosol_absorption_m_inv;
    double ozone_absorption_m_inv;
    double o2_absorption_m_inv;
    double o4_absorption_m_inv;
    double h2o_absorption_m_inv;
    double no2_absorption_m_inv;
    double extinction_m_inv;
    double single_scattering_albedo;
    double aerosol_asymmetry;
};

class WavelengthManager
{
public:
    void configureGrid(double min_wavelength_nm, double max_wavelength_nm, double step_nm);
    void loadSolarSpectrumCsv(const std::string &filename);
    void loadInstrumentResponseCsv(const std::string &filename);
    void loadRayleighCrossSectionCsv(const std::string &filename);
    void loadOzoneCrossSectionCsv(const std::string &filename);
    void loadO2CrossSectionCsv(const std::string &filename);
    void loadO4CrossSectionCsv(const std::string &filename);
    void loadH2OCrossSectionCsv(const std::string &filename);
    void loadNO2CrossSectionCsv(const std::string &filename);
    void buildBands();

    const std::vector<SpectralBand> &getBands() const { return bands_; }
    double solarIrradiance(double wavelength_nm) const;
    double instrumentResponse(double wavelength_nm) const;
    double rayleighCrossSection(double wavelength_nm) const;
    double ozoneCrossSection(double wavelength_nm) const;
    double o2CrossSection(double wavelength_nm) const;
    double o4CrossSection(double wavelength_nm) const;
    double h2oCrossSection(double wavelength_nm) const;
    double no2CrossSection(double wavelength_nm) const;

private:
    double min_wavelength_nm_ = 350.0;
    double max_wavelength_nm_ = 800.0;
    double step_nm_ = 10.0;
    std::vector<std::pair<double, double>> solar_irradiance_table_;
    std::vector<std::pair<double, double>> instrument_response_table_;
    std::vector<std::pair<double, double>> rayleigh_cross_section_table_;
    std::vector<std::pair<double, double>> ozone_cross_section_table_;
    std::vector<std::pair<double, double>> o2_cross_section_table_;
    std::vector<std::pair<double, double>> o4_cross_section_table_;
    std::vector<std::pair<double, double>> h2o_cross_section_table_;
    std::vector<std::pair<double, double>> no2_cross_section_table_;
    std::vector<SpectralBand> bands_;
};

double rayleighCrossSectionM2(double wavelength_nm);
LocalOpticalProperties computeOpticalProperties(
    const Atmosphere &atmosphere,
    const WavelengthManager &wavelength_manager,
    double altitude_m,
    double wavelength_nm
);

#endif // WAVELENGTHHANDLING_HPP
