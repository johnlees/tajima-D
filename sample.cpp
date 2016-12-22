/*
 * File: sample.cpp
 *
 * Helper functions for the sample class
 *
 */

#include "sample.hpp"

Sample::Sample(size_t num_vars)
{
   _seq.zeros(num_vars);
}

// Position in is 1-indexed
void Sample:add_to_seq(double var, size_t pos)
{
   _seq[pos-1] = var;
}

