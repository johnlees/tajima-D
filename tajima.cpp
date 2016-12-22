/*
 * File: tajima.cpp
 *
 * Main control loop for tajima
 *
 */

#include "tajima.hpp"

// Constants
const std::string VERSION = "0.1";
//    Default options

int main (int argc, char *argv[])
{
   // Read cmd line options and process
   // Read in all variant lines, store positions of 1s in vector for each line
   // Create N Sample objects, with var size
   // For each variant line, add any ones into relevant sample at correct
   // positon
   // Using full sample objects, calculate Tajima's D

   // Program description
   std::cerr << "tajima: calculate tajima's D from bcftools query output\n";

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   boost::program_options::variables_map vm;
   if (argc == 1)
   {
      std::cerr << "Usage: epistasis -b bacterial_snps_indels.csv.gz -h human_snps.csv.gz --struct all_structure\n\n"
         << "For full option details run epistasis -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, vm))
   {
      return 1;
   }

   // Open human file to get number of samples
   std::vector<std::string> human_variant;
   if (vm.count("human"))
   {
      igzstream human_in;
      human_in.open(vm["human"].as<std::string>().c_str());
      human_variant = readCsvLine(human_in);
   }
   else
   {
      throw std::runtime_error("--human option is compulsory");
   }

   size_t num_samples = human_variant.size();

   // Read in all the bacterial variants (3Mb compressed - shouldn't be too bad
   // in this form I hope)
   std::cerr << "Reading in all bacterial variants" << std::endl;
   igzstream bacterial_file;
   bacterial_file.open(parameters.bact_file.c_str());

   std::vector<Pair> all_pairs;
   long int bact_line_nr = 1;
   while (bacterial_file)
   {
      std::vector<std::string> bacterial_variant = readCsvLine(bacterial_file);
      if (bacterial_file)
      {
         Pair bact_in(num_samples);
         bact_in.add_y(bacterial_variant, bact_line_nr);

         // Check MAF and missingness of this variant
         std::tuple<double,double> mafs = bact_in.maf();
         std::tuple<double,double> missings = bact_in.missing();
         if (std::get<1>(mafs) > parameters.min_af && std::get<1>(mafs) < parameters.max_af && std::get<1>(missings) < parameters.missing)
         {
            if (use_mds)
            {
               bact_in.add_covar(mds);
            }

            // If set here will calculate logistic regression for all bacterial
            // variants if covar provided.
            // Alternative would be to do for only pairs passing chi-sq. Less
            // efficient if many pairs passing.
            set_null_ll(bact_in);

            all_pairs.push_back(bact_in);
         }

         bact_line_nr++;
      }
      else
      {
         bact_line_nr--;
         break;
      }
   }

   // Write a header
   std::cerr << "Starting association tests" << std::endl;

   while (human_file)
   {
      std::vector<std::string> human_variant;
      human_variant.reserve(num_samples);
      human_variant = readCsvLine(human_file);

      if (human_file)
      {
         // Test each human variant against every bacterial variant
         for (auto it = all_pairs.begin(); it < all_pairs.end(); it++)
         {
            it->add_x(human_variant, human_line_nr);

            // maf filter
            std::tuple<double,double> mafs = it->maf();
            std::tuple<double,double> missings = it->missing();
            if (std::get<0>(mafs) > parameters.min_af && std::get<0>(mafs) < parameters.max_af && std::get<0>(missings) < parameters.missing)
            {
               it->chisq_p(chiTest(*it));
               read_pairs++;

               if (it->chisq_p() < parameters.chi_cutoff)
               {
                  doLogit(*it);

                  // Likelihood ratio test
                  it->p_val(likelihoodRatioTest(*it));

                  tested_pairs++;
                  if (it->p_val() < parameters.log_cutoff)
                  {
                     significant_pairs++;
                  }
               }

               out_stream << *it << std::endl;
            }
         }

         if (parameters.chunk_end > 1 && human_line_nr >= parameters.chunk_end)
         {
            human_line_nr--;
            break;
         }
         else
         {
            human_line_nr++;
         }
      }
      else
      {
         break;
      }
   }

   std::cerr << "Done.\n";
}

std::vector<std::string> readCsvLine(std::istream& is)
{
   std::vector<std::string> variant;
   std::string line;
   std::getline(is, line);

   std::stringstream line_stream(line);
   std::string value;

   while(std::getline(line_stream, value, ','))
   {
      variant.push_back(value);
   }
   return variant;
}

