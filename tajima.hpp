/*
 *
 * tajima.hpp
 * Header file for tajima
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <random>
#include <vector>
#include <exception>
#include <sys/stat.h>

// Boost headers
#include <boost/program_options.hpp>

// Armadillo headers
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

// Constants
extern const std::string VERSION;

// Function headers for each cpp file

// tajima.cpp
std::tuple<long int,std::vector<int>> readCsvLine(std::string& line);
double D_subsample(const arma::mat& alignment, const arma::uvec& sampled_indices, const std::vector<long int>& positions, const int verbose);
double calc_D(const arma::mat& alignment, std::vector<long int>& positions, const int verbose);

// cmdLine.cpp
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);
int fileStat(const std::string& filename);

