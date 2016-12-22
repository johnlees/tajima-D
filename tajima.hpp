/*
 *
 * tajima.hpp
 * Header file for tajima
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <exception>
#include <sys/stat.h>

// Boost headers
#include <boost/program_options.hpp>

// Armadillo headers
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

// Classes
#include "sample.hpp"

// Constants
extern const std::string VERSION;

// Function headers for each cpp file

// tajima.cpp
std::vector<int> readCsvLine(std::istream& is);

// cmdLine.cpp
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);
int fileStat(const std::string& filename);

