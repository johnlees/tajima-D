/*
 * sample.hpp
 * Header file for Sample class
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <iterator>
#include <vector>
#include <exception>

// Armadillo/dlib headers
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

class Sample
{
   public:
      // Initialisation
      Sample(size_t num_vars) {  _seq.zeros(num_vars); };

      // nonmodifying operations
      long int full_seq() const { return _seq; }

      // Modifying operations
      void set_seq(arma::vec full_seq) { _seq = full_seq; }
      void add_to_seq(double var, size_t pos) { _seq[pos] = var; } // defined in sample.cpp

   private:
      arma::vec _seq;

};

