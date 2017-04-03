/*
 * File: cmdLine.cpp
 *
 * Reads command line input to tajima using boost
 * program options
 *
 */

#include "tajima.hpp"

namespace po = boost::program_options; // Save some typing

// Use boost::program_options to parse command line input
// This does pretty much all the parameter checking needed
int parseCommandLine (int argc, char *argv[], po::variables_map& vm)
{
   int failed = 0;

   //Required options
   po::options_description required("Required options");
   required.add_options()
    ("snps,s", po::value<std::string>()->required(), "input snps from bcftools query")

   po::options_description other("Other options");
   other.add_options()
    ("seg_sites", po::value<double>(), "number of segregating sites, from bcftools view -H | wc -l")
    ("verbose,v", "display progress")
    ("version", "prints version and exits")
    ("help,h", "full help message");

   po::options_description all;
   all.add(required).add(other);

   try
   {
      po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

      if (vm.count("help"))
      {
         printHelp(all);
         failed = 1;
      }
      else if (vm.count("version"))
      {
         std::cout << VERSION << std::endl;
         failed = 1;
      }
      else
      {
         po::notify(vm);
         failed = 0;

         // Check input files exist, and can stat
         if (!fileStat(vm["snps"].as<std::string>()))
         {
            failed = 1;
         }
      }

   }
   catch (po::error& e)
   {
      // Report errors from boost library
      std::cerr << "Error in command line input: " << e.what() << "\n";
      std::cerr << "Run 'tajima --help' for full option listing\n\n";
      std::cerr << required << "\n" << other << "\n";

      failed = 1;
   }

   return failed;
}

// Print long help message
void printHelp(po::options_description& help)
{
   std::cerr << help << "\n";
}

// Check for file existence
int fileStat(const std::string& filename)
{
   struct stat buffer;
   int success = 1;

   if (stat (filename.c_str(), &buffer) != 0)
   {
      std::cerr << "Can't stat input file: " << filename << "\n";

      success = 0;
   }

   return success;
}

