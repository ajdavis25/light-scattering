#ifndef VALIDATIONQA_HPP
#define VALIDATIONQA_HPP

#include <string>
#include <vector>

struct ValidationMetric
{
    std::string name;
    double value = 0.0;
    double threshold = 0.0;
    bool pass = true;
};

struct ValidationReport
{
    bool overall_pass = true;
    std::vector<ValidationMetric> metrics;
    std::vector<std::string> notes;
};

ValidationReport runValidationSuite(const std::string &config_path);
void writeValidationReport(const ValidationReport &report, const std::string &output_dir);

#endif // VALIDATIONQA_HPP
