// AccelerationTechniques.hpp
#ifndef ACCELERATIONTECHNIQUES_HPP
#define ACCELERATIONTECHNIQUES_HPP

/**
 * Russian roulette approach:
 * If photon weight < threshold, randomly kill or boost
 */
double russianRoulette(double &photonWeight, double threshold);

/**
 * Photon splitting if it enters high scattering region, etc.
 */
void photonSplitting(/*...*/);

#endif
