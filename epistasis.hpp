/*
 *
 * seercommon.hpp
 * Header file for seercommon
 * Shared functions between seer and kmds
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <chrono>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <thread>
#include <exception>
#include <sys/stat.h>
#include <regex>

// gzstream headers
#include <gzstream.h>

// Boost headers
#include <boost/program_options.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

// Armadillo/dlib headers
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include <dlib/matrix.h>

// dlib headers
#include <dlib/optimization.h>

// Classes
#include "pair.hpp"

// Constants
extern const std::string VERSION;
extern const std::string maf_default;
extern const std::string missing_default;
extern const std::string chisq_default;
extern const std::string pval_default;
extern const double convergence_limit;
extern const unsigned int max_nr_iterations;
extern const double se_limit;
extern const double bfgs_start_beta;

typedef dlib::matrix<double,0,1> column_vector;

// Structs
struct cmdOptions
{
   double log_cutoff;
   double chi_cutoff;

   double min_af;
   double max_af;
   double missing;

   long int chunk_start;
   long int chunk_end;

   std::string bact_file;
   std::string human_file;
   std::string output_file;
   std::string struct_file;
};

// Function headers for each cpp file

// epistasis.cpp
std::vector<std::string> readCsvLine(std::istream& is);

// common.cpp
cmdOptions verifyCommandLine(boost::program_options::variables_map& vm, double num_samples);
arma::vec dlib_to_arma(const column_vector& dlib_vec);
column_vector arma_to_dlib(const arma::vec& arma_vec);
arma::mat inv_covar(arma::mat A);
int fileStat(const std::string& filename);

// cmdLine.cpp
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);

// logisticRegression.cpp
void doLogit(Pair& p);
void newtonRaphson(Pair& p, const arma::vec& y_train, const arma::mat& x_design, const bool firth);
arma::mat varCovarMat(const arma::mat& x, const arma::mat& b);
arma::vec predictLogitProbs(const arma::mat& x, const arma::vec& b);

// stats.cpp
double chiTest(Pair& p);
void set_null_ll(Pair& p);
double likelihoodRatioTest(Pair& p);
double normalPval(double testStatistic);

// fisher.cpp
double fisher22(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t midp);
double fisher22_1sided(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t m11_is_greater_alt, uint32_t midp);
void fisher22_precomp_thresh(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t* m11_minp, uint32_t* m11_maxp, uint32_t* tiep);
int32_t fisher23_tailsum(double* base_probp, double* saved12p, double* saved13p, double* saved22p, double* saved23p, double *totalp, uint32_t* tie_ctp, uint32_t right_side);
double fisher23(uint32_t m11, uint32_t m12, uint32_t m13, uint32_t m21, uint32_t m22, uint32_t m23, uint32_t midp);

