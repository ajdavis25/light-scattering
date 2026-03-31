#include "MonteCarloDriver.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
    try {
        const std::string configPath = argc > 1
            ? argv[1]
            : (std::filesystem::exists("../config/default_clear_sky.cfg")
                ? "../config/default_clear_sky.cfg"
                : "monte_carlo_cpp/config/default_clear_sky.cfg");

        std::cout << "Running production clear-sky twilight solver with config: "
                  << configPath << "\n";
        const SimulationConfig config = loadSimulationConfig(configPath);
        const SkyResult result = runMonteCarloSimulation(config);
        writeSkyResult(result);
        std::cout << "Simulation complete.\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Simulation failed: " << error.what() << "\n";
        return 1;
    }
}
