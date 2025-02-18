// Parallelization.hpp
#ifndef PARALLELIZATION_HPP
#define PARALLELIZATION_HPP

/**
 * Initialize MPI or manage worker distribution
 */
void initParallel();

/**
 * Gather results from each rank
 */
void gatherAndCombine(/*some data structure*/);

#endif
