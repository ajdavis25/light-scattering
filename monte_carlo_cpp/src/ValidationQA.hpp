// ValidationQA.hpp
#ifndef VALIDATIONQA_HPP
#define VALIDATIONQA_HPP

/**
 * Compare single-scatter results with analytic formula for Rayleigh 
 * or compare plane-parallel approach with a known code or reference.
 */
void testSingleScatterPlaneParallel();

/**
 * Test or compare multi-scattering results with known approximate solutions
 * or a simpler reference code.
 */
void testMultiScatterAgainstDISORT();

#endif
