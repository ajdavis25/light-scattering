#include "ValidationQA.hpp"

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
        const ValidationReport report = runValidationSuite(configPath);
        writeValidationReport(report, "monte_carlo_cpp/results/validation");

        std::cout << "Validation overall pass: " << (report.overall_pass ? "true" : "false") << "\n";
        for (const auto &metric : report.metrics) {
            std::cout << metric.name << " = " << metric.value
                      << " (threshold " << metric.threshold << ", pass=" << (metric.pass ? "true" : "false") << ")\n";
        }
        for (const auto &note : report.notes) {
            std::cout << "Note: " << note << "\n";
        }
        return report.overall_pass ? 0 : 1;
    } catch (const std::exception &error) {
        std::cerr << "Validation failed to run: " << error.what() << "\n";
        return 1;
    }
}
